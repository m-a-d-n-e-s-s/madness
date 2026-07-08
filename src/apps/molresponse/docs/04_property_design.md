# Property Design

> **Superseded by `14_property_architecture.md`** for the overall framing.
> The two-tier split (Response Properties = core driver; Chemical Properties
> = post-processing) and the user-facing specification model live there. This
> doc retains only the persistence/matching notes that feed it.

Status: **stub.** Property assembly (polarizability, hyperpolarizability,
resonant Raman, two-photon absorption) is not yet implemented in v3.

The persistence/matching foundation properties depend on is specified in
`13_unified_persistence_schema.md`:

- Properties combining FD and ES inputs (resonant Raman, TPA) are valid only
  when both were computed at the **same protocol** — enforced via the shared
  `protocol_key` and recorded under `properties/<name>/<protocol_key>/` in the
  unified `response_metadata.json`.
- ES-derived frequencies (e.g. the excitation energy used as the FD frequency
  for resonant Raman) are read from `excited_states/<protocol_key>/roots[*].omega`
  for the matching `root_id`, never hard-coded — ω is a per-protocol output.

Fill this in when the first ES-derived property lands (increment 13e).
