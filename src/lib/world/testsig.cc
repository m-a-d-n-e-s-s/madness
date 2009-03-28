#include <iostream>
#include <signal.h>
#include <cstdlib>
#include <cstdio>
#include <execinfo.h>

void print_backtrace(int fd) {
    const int naddr=1024;
    void* buf[naddr];
    int nptrs = backtrace(buf, naddr);
    printf("no. of address from backtrace %d\n", nptrs);
    backtrace_symbols_fd(buf, nptrs, fd);

    // Fox x86-32 do not seem able to get a backtrace from the signal
    // handler even with these compilation options g++ -Wall -g -O0
    // -fno-inline -fasynchronous-unwind-tables
    // -fno-omit-frame-pointer testsig.cc.
    //
    // Have not tried 32-bit.
}

void madness_signal_action(int signum, siginfo_t *info, void *context) {
    if (signum == SIGSEGV) {
        fprintf(stderr,"SIGSEGV: Invalid memory reference ... address=%p\n",info->si_addr);
    }
    else if (signum == SIGILL) {
        fprintf(stderr,"SIGILL: Illegal Instruction ... address=%p\n",info->si_addr);
    }
    else if (signum == SIGFPE) {
        fprintf(stderr,"SIGFPE: Floating point exception ... address=%p\n",info->si_addr);
    }
    else if (signum == SIGBUS) {
        fprintf(stderr,"SIGBUS: Bus error (bad memory access) ... address=%p\n",info->si_addr);
    }
    else {
        fprintf(stderr,"in generic signal handler ... signum=%d\n",signum);
    }

    print_backtrace(fileno(stderr));
    exit(1);
}

void e() {
    void (*f)()=0;
    f();
}
void d() {
    e();
}
void c() {
    d();
}
void b() {
    c();
}
void a() {
    b();
}

int main() {
    struct sigaction act, oldact;
    act.sa_sigaction = madness_signal_action;
    act.sa_flags = SA_SIGINFO;
    if (sigaction(SIGSEGV, &act, &oldact)) printf("failed to install SIGSEGV handler\n");
    if (sigaction(SIGFPE, &act, &oldact)) printf("failed to install SIGSEGV handler\n");
    if (sigaction(SIGILL, &act, &oldact)) printf("failed to install SIGSEGV handler\n");
    if (sigaction(SIGBUS, &act, &oldact)) printf("failed to install SIGSEGV handler\n");

    a();

    return 0;
}
