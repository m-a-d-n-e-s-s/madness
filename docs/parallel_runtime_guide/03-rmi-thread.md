# Chapter 3 — The RMI thread

[← Active Messages](02-active-messages.md) · [Index](README.md) · [Next: TaskQueue & Futures →](04-taskq-futures.md)

Every rank runs exactly **one** RMI (Remote Method Invocation) server thread. It
owns the wire: it posts all receives, polls MPI, dispatches incoming active
messages to their handlers, drains outbound sends queued by other threads, and
runs the huge-message rendezvous. Worker threads never touch MPI for AM — they
hand work to this thread.

Files: `worldrmi.h`, `worldrmi.cc`.

---

## 3.1 Setup: posted receives

At startup (`worldrmi.cc:227-322`) the thread allocates a ring of receive buffers
and posts a non-blocking `Irecv` for each, all on `RMI_TAG` (1023) with
`MPI_ANY_SOURCE`:

- `nrecv_` buffers — `MAD_RECV_BUFFERS`, default 128, min 32.
- each `max_msg_len_` bytes — `MAD_BUFFER_SIZE`, default `3·512·1024` = 1.5 MB,
  min 1024.
- one extra slot (`recv_buf[nrecv_]`) reserved for the current huge message.
- per-source ordering counters `send_counters[]` / `recv_counters[]`.

```
recv_req[0..nrecv_-1]  ── Irecv(buf_i, max_msg_len_, ANY_SOURCE, RMI_TAG)
recv_req[nrecv_]       ── reserved for one in-flight huge message
```

---

## 3.2 The server loop: `process_some()`

This is the heart of the runtime's communication (`worldrmi.cc:61-178`). One pass:

```mermaid
flowchart TB
    start([process_some]) --> poll
    poll["Testsome(recv_req): how many arrived?"] --> any{narrived > 0?}
    any -- no --> backoff["clear_send_req();<br/>usleep(MAD_BACKOFF_US);<br/>retry up to 1000x"]
    backoff --> poll
    any -- yes --> loop[for each arrived message m]
    loop --> read["read header: func, attr;<br/>count = attr>>16;<br/>stats += len"]
    read --> ord{ordered AND<br/>count != recv_counters[src]?}
    ord -- "no (in-order/unordered)" --> run["if ordered: recv_counters[src]++;<br/>func(buf, len);  repost Irecv"]
    ord -- "yes (early arrival)" --> q["push to out-of-order queue q[]"]
    run --> done1
    q --> done1[next message]
    done1 --> drain["sort q by (src,count);<br/>run any now-in-order msgs;<br/>keep leftovers"]
    drain --> huge["post_pending_huge_msg()"]
    huge --> flush["ThreadPool::flush_prebuf();<br/>clear_send_req()"]
    flush --> ret([return])
```

Key points:

1. **Polling, not waiting.** Because MPI is only `THREAD_SERIALIZED`, the thread
   cannot `Waitsome`; it `Testsome`es up to 1000 times, backing off
   `MAD_BACKOFF_US` (default **2 µs**, `worldrmi.cc:499`) between empty polls
   (`worldrmi.cc:76-82`). Between polls it also drains outbound sends with
   `clear_send_req()`.
2. **Dispatch.** For each arrived message it reads the relative function pointer
   and the attribute word, recovers the handler, and — if the message is unordered
   or already in order — calls it immediately, then **reposts** that receive buffer
   (`post_recv_buf`, `worldrmi.cc:201-213`).
3. **Reordering.** An ordered message that arrives early (its `count` ≠ the
   expected `recv_counters[src]`) is parked in an out-of-order queue `q[]`. After
   the arrival loop, `q` is sorted by `(src, count)` — the comparator handles
   16-bit counter wraparound — and any messages that are now next-in-line are run
   (`worldrmi.cc:137-166`). Leftovers stay for the next pass.
4. **Prebuffer flush.** With the Pthreads backend the RMI thread never runs tasks
   or calls `await`, so it must manually flush the thread-local task-submission
   prebuffer (`ThreadPool::flush_prebuf()`, `worldrmi.cc:174`) — otherwise tasks
   that handlers submitted would sit invisibly in a thread-local buffer.

---

## 3.3 The send side: `RMI::isend`

Outbound AMs (queued by `WorldAmInterface::send`) are actually injected here
(`worldrmi.cc:~420-486`):

```cpp
lock();
if (is_ordered(attr)) attr |= ((send_counters[dest]++)<<16);  // stamp sequence
h->func = to_rel_fn_ptr(func);  h->attr = attr;               // write header
stats.nmsg_sent++;  stats.nbyte_sent += nbyte;
numsent++;
if (nssend_ && numsent==nssend_) {            // periodic synchronous send
    result = comm.Issend(buf, nbyte, MPI_BYTE, dest, tag);   // back-pressure
    numsent %= nssend_;
} else {
    result = comm.Isend (buf, nbyte, MPI_BYTE, dest, tag);
}
unlock();
```

The **`MAD_NSSEND` throttle** (default = `nrecv_`) replaces every Nth `Isend` with
a synchronous `Issend` (`worldrmi.cc:475-481`). `Issend` does not complete until
the matching receive is posted, which periodically forces the sender to wait for
the receiver to keep up — a cheap, self-clocking back-pressure mechanism that
prevents a fast sender from overrunning a slow receiver's buffer ring.

---

## 3.4 Huge messages (rendezvous)

A message larger than `max_msg_len_` cannot land in a normal receive buffer, so it
uses a rendezvous on the 4096–8191 tag range:

```mermaid
sequenceDiagram
    autonumber
    participant S as sender (isend)
    participant R as receiver RMI thread
    S->>S: tag = unique_tag()  (4096..8191)
    S->>R: control AM (huge_msg_handler): {src, nbyte, tag}  [UNORDERED]
    S->>S: Irecv(ack, src, tag + unique_tag_period())
    R->>R: queue (src,nbyte,tag) in hugeq
    Note over R: when the reserved huge slot is free
    R->>R: posix_memalign(buf, nbyte); Irecv(buf, nbyte, src, tag)
    R->>S: Send(ack, src, tag + unique_tag_period())
    S->>R: Isend(payload, nbyte, dest, tag)
    Note over R: payload lands in the huge buffer; handler runs; buffer freed
```

This is `post_pending_huge_msg` (`worldrmi.cc:180-199`) plus `huge_msg_handler`.
Only one huge message is in flight per rank at a time (the single reserved slot),
so frequent huge messages serialize. **Raising `MAD_BUFFER_SIZE`** so that typical
coefficient tensors fit inline avoids this path entirely (Chapter 9).

---

## 3.5 Statistics you can read

`RMI::stats` (`worldrmi.cc:51, 96-98, 469-470`) tracks `nmsg_sent`, `nbyte_sent`,
`nmsg_recv`, `nbyte_recv`, and `max_serv_send_q`. These are the raw inputs to the
communication term of any performance model (companion doc §9). `World::print_stats`
(`world.cc:236-421`) dumps them at shutdown.

---

## 3.6 Latency anatomy (for modeling)

A single AM's wall-clock latency has three additive parts:

```
 t_AM ≈ t_send_buffer_wait      (0 unless > nsend in flight → 100µs sleeps)
      + t_mpi_transit           (network; ~Issend stall every MAD_NSSEND msgs)
      + t_receiver_poll+handler (≤ one MAD_BACKOFF_US interval + handler time)
```

- Latency-bound, many tiny messages → raise `MAD_SEND_BUFFERS`, lower
  `MAD_BACKOFF_US`.
- Bandwidth-bound, large tensors → raise `MAD_BUFFER_SIZE` (avoid rendezvous), and
  remember the wire is single-threaded (Chapter 1).

[← Active Messages](02-active-messages.md) · [Index](README.md) · [Next: TaskQueue & Futures →](04-taskq-futures.md)
