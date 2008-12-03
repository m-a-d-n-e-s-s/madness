#ifndef MAD_DISPL_H
#define MAD_DISPL_H

namespace madness {
    /// Holds displacements for applying operators to avoid replicating for all operators
    template <int NDIM>
    class Displacements {

        static std::vector< Key<NDIM> > disp;
        static std::vector< Key<NDIM> > disp_periodicsum[64];

        static int bmax_default() {
            int bmax;
            if (NDIM == 1) bmax = 7;

            else if (NDIM == 3) bmax = 3;
            else bmax = 2;
            return bmax;
        }

        static bool cmp_keys(const Key<NDIM>& a, const Key<NDIM>& b) {
            return a.distsq() < b.distsq();
        }

        static bool cmp_keys_periodicsum(const Key<NDIM>& a, const Key<NDIM>& b) {
            Translation twonm1 = (Translation(1)<<a.level())>>1;

            uint64_t suma=0, sumb=0;
            for (int d=0; d<NDIM; d++) {
                Translation la = a.translation()[d];
                if (la > twonm1) la -= twonm1*2;
                if (la <-twonm1) la += twonm1*2;
                suma += la*la;

                Translation lb = b.translation()[d];
                if (lb > twonm1) lb -= twonm1*2;
                if (lb <-twonm1) lb += twonm1*2;
                sumb += lb*lb;
            }
            return suma < sumb;
        }

        static void make_disp(int bmax) {
            // Note newer loop structure in make_disp_periodic_sum
            Vector<Translation,NDIM> d;

            int num = 1;
            for (int i=0; i<NDIM; i++) num *= (2*bmax + 1);
            disp = std::vector< Key<NDIM> >(num);

            num = 0;
            if (NDIM == 1) {
                for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                    disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 2) {
                for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                    for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                        disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 3) {
                for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                    for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                        for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                            disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 4) {
                for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                    for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                        for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                            for (d[3]=-bmax; d[3]<=bmax; d[3]++)
                                disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 5) {
                for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                    for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                        for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                            for (d[3]=-bmax; d[3]<=bmax; d[3]++)
                                for (d[4]=-bmax; d[4]<=bmax; d[4]++)

                                    disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 6) {
                for (d[0]=-bmax; d[0]<=bmax; d[0]++)
                    for (d[1]=-bmax; d[1]<=bmax; d[1]++)
                        for (d[2]=-bmax; d[2]<=bmax; d[2]++)
                            for (d[3]=-bmax; d[3]<=bmax; d[3]++)
                                for (d[4]=-bmax; d[4]<=bmax; d[4]++)
                                    for (d[5]=-bmax; d[5]<=bmax; d[5]++)
                                        disp[num++] = Key<NDIM>(0,d);
            }
            else {
                MADNESS_EXCEPTION("_make_disp: hard dimension loop",NDIM);
            }

            std::sort(disp.begin(), disp.end(), cmp_keys);
        }

        static void make_disp_periodicsum(int bmax, Level n) {
            Translation twon = Translation(1)<<n;

            if (bmax > (twon-1)) bmax=twon-1;

            // Make permissible 1D translations
            Translation b[4*bmax+1];
            int i=0;
            for (Translation lx=-bmax; lx<=bmax; lx++) {
                b[i++] = lx;
                if ((lx < 0) && (lx+twon > bmax)) b[i++] = lx + twon;
                if ((lx > 0) && (lx-twon <-bmax)) b[i++] = lx - twon;
            }
            MADNESS_ASSERT(i <= 4*bmax+1);
            int numb = i;

            disp_periodicsum[n] = std::vector< Key<NDIM> >();
            Vector<long,NDIM> lim(numb);
            for (IndexIterator index(lim); index; ++index) {
                Vector<Translation,NDIM> d;
                for (int i=0; i<NDIM; i++) {
                    d[i] = b[index[i]];
                }
                disp_periodicsum[n].push_back(Key<NDIM>(n,d));
            }

            std::sort(disp_periodicsum[n].begin(), disp_periodicsum[n].end(), cmp_keys_periodicsum);
//             print("KEYS AT LEVEL", n);
//             print(disp_periodicsum[n]);
        }


    public:
        Displacements() {
            if (disp.size() == 0) {
                make_disp(bmax_default());

                Level nmax = 8*sizeof(Translation) - 2;
                for (Level n=0; n<nmax; n++) make_disp_periodicsum(bmax_default(), n);
            }
        }

        const std::vector< Key<NDIM> >& get_disp(Level n, bool isperiodicsum) {
            if (isperiodicsum) {
                return disp_periodicsum[n];
            }
            else {
                return disp;
            }
        }

    };
}
#endif
