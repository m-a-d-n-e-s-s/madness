#ifndef MAD_PRINT_SEQ_H
#define MAD_PRINT_SEQ_H

/// \file print_seq.h
/// \brief Implements print_seq ... included by world.h

namespace madness {
    /// Sequentially ordered printing of (serializable) data from every process ... collective no fence
    template <typename A, typename B, typename C>
    void print_seq(World& world, const A& a, const B& b, const C& c) {
        if (world.rank() == 0) {
            printf("%6d : ",0);
            print(a, b, c);
            for (int p=1; p<world.size(); p++) {
                A aa;
                B bb;
                C cc;
                MPIOutputArchive(world,p) & 1;
                MPIInputArchive(world, p) & aa & bb & cc;
                printf("%6d : ",p);
                print(aa, bb, cc);
            }
        }
        else {
            int i;
            MPIInputArchive(world, 0) & i;
            MPIOutputArchive(world, 0) & a & b & c;
        }
    }

    /// Sequentially ordered printing of (serializable) data from every process ... collective no fence
    template <typename A, typename B>
    void print_seq(World& world, const A& a, const B& b) {
        if (world.rank() == 0) {
            printf("%6d : ",0);
            print(a, b);
            for (int p=1; p<world.size(); p++) {
                A aa;
                B bb;
                MPIOutputArchive(world,p) & 1;
                MPIInputArchive(world, p) & aa & bb;
                printf("%6d : ",p);
                print(aa, bb);
            }
        }
        else {
            int i;
            MPIInputArchive(world, 0) & i;
            MPIOutputArchive(world, 0) & a & b;
        }
    }

    /// Sequentially ordered printing of (serializable) data from every process ... collective no fence
    template <typename A>
    void print_seq(World& world, const A& a) {
        if (world.rank() == 0) {
            printf("%6d : ",0);
            print(a);
            for (int p=1; p<world.size(); p++) {
                A aa;
                MPIOutputArchive(world,p) & 1;
                MPIInputArchive(world, p) & aa;
                printf("%6d : ",p);
                print(aa);
            }
        }
        else {
            int i;
            MPIInputArchive(world, 0) & i;
            MPIOutputArchive(world, 0) & a;
        }
    }
}

#endif
