""" Codegen for mtxm

A   I
  +-------------+
K | a b ....... |
  |     ...     |
  +-------------+

B   J
  +-------------+
K | i j k l ... |
  |     ...     |
  +-------------+

C   J
  +-------------+
I | w x y z ... |
  |     ...     |
  +-------------+

# Complex Complex
w += a i - b j
x += a j + b i
y += a k - b l
z += a l + b k

temps: _c += _a _b +- _ai _br


# Real Complex
w += a i
x += a j
y += a k
z += a l

temps: _c += _a _b ; doubled dimj


# Complex Real
w += a i
x += b i
y += a j
z += b j

temps: _c += _az _bz


# Real Real
w += a i
x += a j
y += a k
z += a l

temps: _c += _a _b
"""

from itertools import product
import logging
logger = logging.getLogger(__name__)

class MTXMGen:
    def __init__(self, cxa=False, cxb=False):
        self.indent = 4
        self.__in_main_loop = False
        self._mask = False
        self._odds = [1]
        self.have_bgp = False
        self.have_bgq = False
        self.have_avx2 = False
        self.complex_a = cxa
        self.complex_b = cxb
        self.complex_dup_cast = ''

    @property
    def complex_c(self):
        return self.complex_a or self.complex_b

    @property
    def complex_complex(self):
        return self.complex_a and self.complex_b

    @property
    def real_complex(self):
        return not self.complex_a and self.complex_b

    @property
    def complex_real(self):
        return self.complex_a and not self.complex_b

    @property
    def real_real(self):
        return not self.complex_a and not self.complex_b

    def _temp(self, prefix, x, y):
        return prefix + '_' + str(x) + '_' + str(y)

    def _temps(self, prefix, x, y, size):
        """
        >>> list(MTXMGen()._temps('_x', 'i', 'j', {'i':2, 'j':3}))
        ['_x_0_0', '_x_0_1', '_x_0_2', '_x_1_0', '_x_1_1', '_x_1_2']
        """
        return [self._temp(prefix, i, j) for i, j in product(range(size[x]), range(size[y]))]

    def _post_process(self, lines):
        return lines

    def _header(self, func_name):
        f = lambda x: x and "complex" or ""
        x2 = lambda x: x and "*2" or ""
        if self.have_bgq:
            ret = ["void " + func_name + "(long dimi, long dimj, long dimk, long extb, double {} * __restrict__ c_x, const double {} * __restrict__ a_x, const double {} * __restrict__ b_x) {{".format(f(self.complex_c), f(self.complex_a), f(self.complex_b))]
        else:
            ret = ["void " + func_name + "(long dimi, long dimj, long dimk, double {} * __restrict__ c_x, const double {} * __restrict__ a_x, const double {} * __restrict__ b_x) {{".format(f(self.complex_c), f(self.complex_a), f(self.complex_b))]
        ret.append("    int i, j, k;")
        ret.append("    double * __restrict__ c = (double*)c_x;")
        ret.append("    const double * __restrict__ a = (double*)a_x;")
        ret.append("    const double * __restrict__ b = (double*)b_x;")
        ret.append("    long effj = dimj;")
        if self.have_bgq:
            ret.append("    double * c_buf;")
            ret.append("    double * b_buf;")
            ret.append("    bool free_b = false;")
            ret.append("    /* Setup a buffer for c if needed */")
            ret.append("    double* c_out = c;")
            ret.append("    if (dimj%{}) {{".format(self.complex_c and 2 or 4))
            ret.append("        effj = (dimj | {}) + 1;".format(self.complex_c and 1 or 3))
            ret.append("        posix_memalign((void **) &c, 32, dimi*effj*sizeof(double){});".format(x2(self.complex_c)))
            ret.append("        c_buf = c;")
            ret.append("    }")
            ret.append("    /* Copy b into a buffer if needed */")
            ret.append("    if (extb%{}) {{".format(self.complex_b and 2 or 4))
            ret.append("        long t_extb = (dimj | {}) + 1;".format(self.complex_b and 1 or 3))
            ret.append("        free_b = true;")
            ret.append("        posix_memalign((void **) &b_buf, 32, dimk*t_extb*sizeof(double){});".format(x2(self.complex_b)))
            ret.append("        double* bp = b_buf;")
            ret.append("        for (k=0; k<dimk; k++, bp += t_extb{0}, b += extb{0})".format(x2(self.complex_b)))
            ret.append("            memcpy(bp, b, sizeof(double)*dimj{});".format(x2(self.complex_b)))
            ret.append("        b = b_buf;")
            ret.append("        extb = t_extb;")
            ret.append("    }")
        return ret

    def _footer(self):
        x2 = lambda x: x and "*2" or ""
        ret = []
        if self.have_bgq:
            ret.append("    /* Copy c out if needed */")
            ret.append("    if (dimj%{}) {{".format(self.complex_c and 2 or 4))
            ret.append("        double* ct = c_buf;")
            ret.append("        for (i=0; i<dimi; i++, ct += effj{0}, c_out += dimj{0})".format(x2(self.complex_c)))
            ret.append("            memcpy(c_out, ct, sizeof(double)*dimj{});".format(x2(self.complex_c)))
            ret.append("        free(c_buf);")
            ret.append("    }")
            ret.append("    /* Free the buffer for b */")
            ret.append("    if (free_b) free(b_buf);")
        ret.append("}")
        return ret

    def _temp_dec(self, size):
        ret = []
        indent = '    '
        x = indent + self.vector_type + ' '
        x += ', '.join(self._temps('_c', 'i', 'j', size) +
                self._temps('_b', 'k', 'j', size)) + ';'
        ret.append(x)
        x = indent + self.splat_type + ' ' + ', '.join(self._temps('_a', 'k', 'i', size)) + ';'
        ret.append(x)
        if self.complex_complex:
            if not self.have_bgp and not self.have_bgq:
                # BGP does not need seperate reversed registers because a special fma is used
                x = "{} {} {};".format(indent, self.vector_type, ', '.join(self._temps('_br', 'k', 'j', size)))
                ret.append(x)
            if not self.have_bgq:
                # Imaginary component of A
                x = "{} {} {};".format(indent, self.splat_type, ', '.join(self._temps('_ai', 'k', 'i', size)))
                ret.append(x)
        elif self.complex_real:
            # register from A: a b a b
            x = "{} {} {};".format(indent, self.vector_type, ', '.join(self._temps('_az', 'k', 'i', size)))
            ret.append(x)
            # register from B: i i j j
            x = "{} {} {};".format(indent, self.splat_type, ', '.join(self._temps('_bz', 'k', 'j', size)))
            ret.append(x)
        elif self.real_complex:
            pass
        return ret

    def _extra(self):
        return []

    def _temps_to_load(self, unrolls, z, x, y, tname=None):
        if not tname:
            tname = '_'+z
        ret = []
        ystep = 1
        if y == 'j':
            ystep = self.vector_length
        for i, j in product(range(unrolls[x]), range(0, unrolls[y], ystep)):
            ret.append((self._temp(tname, i, j), i, j))
        return ret

    def _load_a(self, unrolls, indent):
        spaces = ' ' * (self.indent*indent)
        ret = []
        for temp, k, i in self._temps_to_load(unrolls, 'a', 'k', 'i'):
            addr = '(pa+' + str((self.complex_a and 2 or 1)*i) + ')'
            if self.complex_real:
                ret.append(self._load_az(spaces, addr, temp, k, i))
            elif self.complex_complex and self.have_bgq:
                ret.append(spaces + temp + " = vec_ld2(0, {});".format(addr))
            else:
                arg0 = ''
                if self.have_bgq:
                    arg0 = '0, '
                ret.append(spaces + temp + ' = {}({}{});'.format(self.splat_op, arg0, addr))
                if self.complex_complex:
                    ret.append(spaces + self._temp('_ai', k, i) + ' = {}({}{}+1);'.format(self.splat_op, arg0, addr))
        return ret

    def _load_b(self, unrolls, indent):
        spaces = ' ' * (self.indent*indent)
        ret = []
        for temp, k, j in self._temps_to_load(unrolls, 'b', 'k', 'j'):
            arg0 = ""
            if self.have_bgq:
                arg0 = "0, "
            addr = '({}pb+{})'.format(arg0, j // (self.complex_real and 2 or 1))
            if self.complex_real:
                ret.append(self._load_bz(spaces, addr, temp, k, j))
            else:
                ret.append(spaces + temp + ' = ' + self.vector_load + addr + ';')
                if self.complex_complex and not self.have_bgp and not self.have_bgq:
                    ret.append(self._load_br(spaces, addr, temp, k, j))
        return ret

    def _load_c(self, unrolls, indent):
        spaces = ' ' * (self.indent*indent)
        ret = []
        for temp, i, j in self._temps_to_load(unrolls, 'c', 'i', 'j'):
            ret.append(spaces + temp + ' = ' + self.vector_zero + ';')
        return ret

    def _load_br(self, spaces, addr, temp, k, j):
        return spaces + self._temp('_br', k, j) + ' = {}({});'.format(self.complex_reverse_dup, addr)

    def _load_az(self, spaces, addr, temp, k, i):
        arg0 = ''
        if self.have_bgq:
            arg0 = '0, '
        return spaces + self._temp('_az', k, i) + ' = {}({}{}{});'.format(self.complex_dup, arg0, self.complex_dup_cast, addr)

    def _load_bz(self, spaces, addr, temp, k, j):
        return spaces + self._temp('_bz', k, j) + ' = {}({});'.format(self.pair_splat, addr)

    def _fma(self, at, bt, ct):
        raise NotImplementedError()

    def _fmaddsub(self, at, bt, ct):
        raise NotImplementedError()

    def _maths(self, unrolls, indent=0):
        spaces = ' ' * (self.indent*indent)
        ret = []
        for j, i, k in product(range(0, unrolls['j'], self.vector_length), range(unrolls['i']), range(unrolls['k'])):
            if self.real_real or self.real_complex:
                at = self._temp('_a', k, i)
                bt = self._temp('_b', k, j)
                ct = self._temp('_c', i, j)
                ret.append(spaces + self._fma(at, bt, ct))
            elif self.complex_real:
                at = self._temp('_az', k, i)
                bt = self._temp('_bz', k, j)
                ct = self._temp('_c', i, j)
                ret.append(spaces + self._fma(at, bt, ct))
            elif self.complex_complex:
                if self.have_avx2:
                    at = self._temp('_ai', k, i)
                    bt = self._temp('_br', k, j)
                    ct = self._temp('_c', i, j)
                    ret.append(spaces + self._fmaddsub(at, bt, ct))
                    at = self._temp('_a', k, i)
                    bt = self._temp('_b', k, j)
                    ret.append(spaces + self._fmaddsub(at, bt, ct))
                else:
                    at = self._temp('_a', k, i)
                    bt = self._temp('_b', k, j)
                    ct = self._temp('_c', i, j)
                    ret.append(spaces + self._fma(at, bt, ct))
                    if not self.have_bgq:
                        at = self._temp('_ai', k, i)
                        if not self.have_bgp:
                            bt = self._temp('_br', k, j)
                    ret.append(spaces + self._fmaddsub(at, bt, ct))
        return ret

    def _array(self, z, x, xx, y, yy, cpx):
        if y == 'j':
            y = "effj"
        else:
            y = "dim" + y

        return z + '+(' + x + '+' + xx + ')*' + y + (cpx and "*2" or "") + '+' + yy

    def _store_c(self, unrolls, indent, bc_mod=""):
        spaces = ' ' * (self.indent*indent)
        ret = []
        jstep = self.vector_length
        for i, j in product(range(unrolls['i']), range(0, unrolls['j'], jstep)):
            if j + jstep < unrolls['j'] or self.__in_main_loop or not self._mask:
                arg0 = self._array(bc_mod+'c', 'i', str(i), 'j', str(j), self.complex_c)
                arg1 = self._temp('_' + 'c', i, j)
                mid = ', '
                if self.have_bgq:
                    arg0, arg1 = arg1, arg0
                    mid = ', 0, '
                ret.append(spaces + self.vector_store + '(' + arg0 + mid + arg1 + ');')
            else:
                # This is somewhat AVX specific, but no other arch's currently support masking, so ok.
                ret.append(spaces + '{}('.format(self.mask_store) + self._array(bc_mod+'c', 'i', str(i), 'j', str(j), self.complex_c) + ', mask, ' + self._temp('_' + 'c', i, j) + ');')
        return ret

    def _loops(self, i, size, bc_mod=""):
        if i == 'i':
            start = 'i=0'
            #FIXME Don't include _odds if i%2==0 and only evens
            loops = [size[i]]
            if loops[-1] != 1:
                loops += self._odds
            if self.have_bgp:
                loops = range(size[i], 0, -2)
            for loop in loops:
                yield ('for ({0}; i+{1}<=dimi; i+={1}) {{'.format(start, loop), loop)
                start = ''
        elif i == 'j':
            loop = size[i] // (self.complex_c and 2 or 1)
            self.__in_main_loop = True
            yield ("for (j=effj; j>{0}; j-={0},{1}c+={0}{2},{1}b+={0}{3}) {{".format(loop, bc_mod, self.complex_c and "*2" or "", self.complex_b and "*2" or ""), size[i])
            self.__in_main_loop = False
            start = ''
            for loop in range(size[i]-self.vector_length, 0, -self.vector_length):
                yield (start + "if (j>{}) {{".format(loop//(self.complex_c and 2 or 1)), loop+self.vector_length)
                start = 'else '
            if size[i] == self.vector_length:
                yield ("{", self.vector_length)
            else:
                yield ("else {", self.vector_length)
        elif i == 'k':
            assert(size[i] == 1)
            pb_inc  = 'effj'
            if self.have_bgq:
                pb_inc = 'extb'
            yield ("for (k=0; k<dimk; k+=1,pb+={}{},pa+=dimi{}) {{".format(pb_inc, self.complex_b and "*2" or "", self.complex_a and "*2" or ""), 1)

    def _close_braces(self, indent=0):
        ret = []
        for i in range(indent, -1, -1):
            ret += [' '*(self.indent*i) + '}']
        return ret

    def _inner_loops(self, perm, sizes, indent=0, unrolls=None, bc_mod=""):
        indent += 1
        if not unrolls:
            unrolls = {x:0 for x in perm}
        ret = []
        spaces = ' '*(self.indent*indent)

        if perm == ['k']:
            ret.append(spaces + "const double* __restrict__ pb = {}b;".format(bc_mod))
            ret.append(spaces + "const double* __restrict__ pa = a+i{};".format(self.complex_a and "*2" or ""))
            ret += self._load_c(unrolls, indent)

        if perm == ['j', 'k']:
            bc_mod = "x"
            ret.append(spaces + "const double* __restrict__ {}b = b;".format(bc_mod))
            ret.append(spaces + "double* __restrict__ {}c = c;".format(bc_mod))

        for loop, unroll in self._loops(perm[0], sizes, bc_mod):
            unrolls[perm[0]] = unroll
            ret.append(spaces + loop)
            if len(perm) > 1:
                ret += self._inner_loops(perm[1:], sizes, indent, unrolls, bc_mod)
            else:
                ret += self._load_a(unrolls, indent+1)
                b_loads = self._load_b(unrolls, indent+1)
                maths = self._maths(unrolls, indent+1)
                b_take = (self.complex_complex and not self.have_bgp and not self.have_bgq) and 2 or 1
                m_take = unrolls['i']*(self.complex_complex and 2 or 1)
                while b_loads:
                    ret += b_loads[0:b_take]
                    ret += maths[0:m_take]
                    b_loads = b_loads[b_take:]
                    maths = maths[m_take:]
            ret.append(spaces + '}')

        if perm == ['k']:
            ret += self._store_c(unrolls, indent, bc_mod)

        return ret

    def gen(self, f, perm, size, func_name='mtxmq'):
        """Output generated code to file f

        Input:
            perm - an array of 'i', 'j', 'k' in the desired loop order
            size - { index : int, }

        Output:
            None
            Code printed to file f
        """

        if type(perm) is not list:
            perm = list(perm)

        if perm[-1] != 'k':
            raise Exception("k must be inner loop")

        lines = []

        # Header
        lines += self._header(func_name)

        # Temps Declaration
        lines += self._temp_dec(size)

        # Architecture Specific declarations, e.g. mask prep
        lines += self._extra()

        # Computation
        lines += self._inner_loops(perm, size)

        # Footer
        lines += self._footer()

        lines = self._post_process(lines)

        # Output
        for line in lines:
            print(line, file=f)


class MTXMAVX(MTXMGen):
    def __init__(self, *args):
        super().__init__(*args)
        self.vector_length = 4

        self.vector_type = '__m256d'
        self.vector_load = '_mm256_loadu_pd'
        self.vector_store = '_mm256_storeu_pd'
        self.vector_zero = '_mm256_setzero_pd()'

        self._mask = True
        self.mask_store = '_mm256_maskstore_pd'

        self.splat_type = '__m256d'
        self.splat_op = '_mm256_broadcast_sd'

        #self.complex_reverse_dup = '_mm256_permute_pd' # (_mm256_loadu_pd(addr), 5), could also use shuffle 5
        self.complex_dup = '_mm256_broadcast_pd'
        self.complex_dup_cast = '(const __m128d*)'
        #self.pair_splat = ''

    def _load_bz(self, spaces, addr, temp, k, j):
        return spaces + self._temp('_bz', k, j) + ' = _mm256_permute_pd(_mm256_broadcast_pd((const __m128d*){}),12);'.format(addr)

    def _load_br(self, spaces, addr, temp, k, j):
        return spaces + self._temp('_br', k, j) + ' = _mm256_permute_pd({}, 5);'.format(temp)

    def _fma(self, at, bt, ct):
        return ct + ' = _mm256_add_pd(_mm256_mul_pd(' + bt + ', ' + at + '), ' + ct + ');'

    def _fmaddsub(self, at, bt, ct):
        return ct + ' = _mm256_addsub_pd(' + ct + ', _mm256_mul_pd(' + at + ', ' + bt + '));'

    def _extra(self):
        if self.real_real:
            return [' ' * self.indent + """
    __m256i mask;
    j = effj % 4;
    switch (j) {
        case 0:
            mask = _mm256_set_epi32(-1,-1,-1,-1,-1,-1,-1,-1);
            break;
        case 1:
            mask = _mm256_set_epi32( 0, 0, 0, 0, 0, 0,-1,-1);
            break;
        case 2:
            mask = _mm256_set_epi32( 0, 0, 0, 0,-1,-1,-1,-1);
            break;
        case 3:
            mask = _mm256_set_epi32( 0, 0,-1,-1,-1,-1,-1,-1);
            break;
        default:
            return;
    }"""]
        else:
            return [' ' * self.indent + """
    __m256i mask;
    j = effj % 2;
    switch (j) {
        case 0:
            mask = _mm256_set_epi32(-1,-1,-1,-1,-1,-1,-1,-1);
            break;
        case 1:
            mask = _mm256_set_epi32( 0, 0, 0, 0,-1,-1,-1,-1);
            break;
        default:
            return;
    }"""]

class MTXMAVX2(MTXMAVX):
    def __init__(self, *args):
        super().__init__(*args)
        self.have_avx2 = True

    def _fma(self, at, bt, ct):
        return ct + ' = _mm256_fmadd_pd(' + ','.join([at,bt,ct]) + ');'

    def _fmaddsub(self, at, bt, ct):
        return ct + ' = _mm256_fmaddsub_pd(' + ','.join([at,bt,ct]) + ');'

class MTXMSSE(MTXMGen):
    def __init__(self, *args):
        super().__init__(*args)
        self.vector_length = 2

        self.vector_type = '__m128d'
        self.vector_load = '_mm_loadu_pd'
        self.vector_store = '_mm_storeu_pd'
        self.vector_zero = '_mm_setzero_pd()'

        self.splat_type = '__m128d'
        self.splat_op = '_mm_load1_pd'

        self.complex_reverse_dup = '_mm_loadr_pd' # aligned only!
        self.complex_dup = '_mm_loadu_pd'
        self.pair_splat = '_mm_load1_pd'

    def _fma(self, at, bt, ct):
        return ct + ' = _mm_add_pd(_mm_mul_pd(' + bt + ', ' + at + '), ' + ct + ');'

    def _fmaddsub(self, at, bt, ct):
        return "{2} = _mm_addsub_pd({2}, _mm_mul_pd({0}, {1}));".format(at, bt, ct)


class MTXMBGP(MTXMGen):
    def __init__(self, *args):
        super().__init__(*args)
        self.have_bgp = True
        self.vector_length = 2

        self.vector_type = '__complex__ double'
        self.vector_load = '__lfpd'
        self.vector_store = '__stfpd'
        self.vector_zero = '__cmplx(0.0,0.0)'

        self.splat_type = 'double'
        self.splat_op = '*'

        self.complex_reverse_dup = '__lfxd'
        self.complex_dup = '__lfpd'
        self.pair_splat = '*'

    def _fma(self, at, bt, ct):
        return ct + ' = __fxcpmadd(' + ct + ', ' + bt + ', ' + at + ');'

    def _fmaddsub(self, at, bt, ct):
        return ct + ' = __fxcxnpma(' + ct + ', ' + bt + ', ' + at + ');'

    def _post_process(self, lines):
        return [x.replace("__restrict__", "").replace("const", "").replace("double complex", "__complex__ double") for x in lines]

class MTXMBGQ(MTXMGen):
    def __init__(self, *args):
        super().__init__(*args)
        self.have_bgq = True
        self.vector_length = 4

        self.vector_type = 'vector4double'
        self.vector_load = 'vec_ld'
        self.vector_store = 'vec_st'
        self.vector_zero = '(vector4double)(0.0)'

        self.complex_dup = 'vec_ld2'

        self.splat_type = 'vector4double'
        self.splat_op = 'vec_lds'

    def _load_bz(self, spaces, addr, temp, k, j):
        t = self._temp('_bz', k, j)
        ret = spaces + t + ' = vec_ld2{};\n'.format(addr)
        ret += spaces + t + ' = vec_perm({0}, {0}, _cr_perm);'.format(t)
        return ret

    def _fma(self, at, bt, ct):
        if self.complex_complex:
            return ct + ' = vec_xmadd(' + at + ', ' + bt + ', ' + ct + ');'
        else:
            return ct + ' = vec_madd(' + at + ', ' + bt + ', ' + ct + ');'

    def _fmaddsub(self, at, bt, ct):
        return ct + ' = vec_xxnpmadd(' + bt + ', ' + at + ', ' + ct + ');'

    def _post_process(self, lines):
        return [x.replace("__restrict__", "").replace("const", "").replace("double complex", "__complex__ double") for x in lines]

    def _extra(self):
        if self.complex_real:
            return [' ' * self.indent + "vector4double _cr_perm = vec_gpci(0x9);"]
        else:
            return []
