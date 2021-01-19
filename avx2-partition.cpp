#include <x86intrin.h>
#include <immintrin.h>
#include <cstdint>


namespace qs {

    namespace avx2 {

        //****************************************************************
        //32 bitni integer brojevi
        __m256i FORCE_INLINE bitmask_to_bytemask_epi32(uint8_t bm) {

            const __m256i mask = _mm256_set1_epi32(bm);
            const __m256i bits = _mm256_setr_epi32(0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80);
            const __m256i tmp = _mm256_and_si256(mask, bits);

            return _mm256_cmpeq_epi32(tmp, bits);
        }


        void FORCE_INLINE
        align_masks(uint8_t &a, uint8_t &b, uint8_t &rem_a, uint8_t &rem_b, __m256i &shuffle_a, __m256i &shuffle_b) {

            assert(a != 0);
            assert(b != 0);

            uint8_t tmpA = a;
            uint8_t tmpB = b;

            uint32_t __attribute__((__aligned__(32))) tmpshufa[8];
            uint32_t __attribute__((__aligned__(32))) tmpshufb[8];

            while (tmpA != 0 && tmpB != 0) {
                int idx_a = __builtin_ctz(tmpA);
                int idx_b = __builtin_ctz(tmpB);

                tmpA = tmpA & (tmpA - 1);
                tmpB = tmpB & (tmpB - 1);

                tmpshufb[idx_a] = idx_b;
                tmpshufa[idx_b] = idx_a;
            }

            a = a ^ tmpA;
            b = b ^ tmpB;

            assert(a != 0);
            assert(b != 0);
            assert(_mm_popcnt_u64(a) == _mm_popcnt_u64(b));

            rem_a = tmpA;
            rem_b = tmpB;

            shuffle_a = _mm256_load_si256((__m256i *) tmpshufa);
            shuffle_b = _mm256_load_si256((__m256i *) tmpshufb);
        }


        __m256i FORCE_INLINE merge(const __m256i mask, const __m256i a, const __m256i b) {

            return _mm256_or_si256(
                    _mm256_and_si256(mask, a),
                    _mm256_andnot_si256(mask, b)
            );
        }


        void FORCE_INLINE swap_epi32(
                __m256i &a, __m256i &b,
                uint8_t mask_a, const __m256i shuffle_a,
                uint8_t mask_b, const __m256i shuffle_b) {

            const __m256i to_swap_b = _mm256_permutevar8x32_epi32(a, shuffle_a);
            const __m256i to_swap_a = _mm256_permutevar8x32_epi32(b, shuffle_b);
            const __m256i ma = bitmask_to_bytemask_epi32(mask_a);
            const __m256i mb = bitmask_to_bytemask_epi32(mask_b);

            a = merge(ma, to_swap_a, a);
            b = merge(mb, to_swap_b, b);
        }


#define _mm256_iszero(vec) (_mm256_testz_si256(vec, vec) != 0)

        void FORCE_INLINE partition_epi32(uint32_t *array, uint32_t pv, int &left, int &right) {

            const int N = 8; // the number of items in a register (256/32)

            __m256i L;
            __m256i R;
            uint8_t maskL = 0;
            uint8_t maskR = 0;

            const __m256i pivot = _mm256_set1_epi32(pv);

            int origL = left;
            int origR = right;

            while (true) {

                if (maskL == 0) {
                    while (true) {
                        if (right - (left + N) + 1 < 2 * N) {
                            goto end;
                        }

                        L = _mm256_loadu_si256((__m256i *) (array + left));
                        const __m256i bytemask = _mm256_cmpgt_epi32(pivot, L);

                        if (_mm256_testc_ps((__m256) bytemask, (__m256) _mm256_set1_epi32(-1))) {
                            left += N;
                        } else {
                            maskL = ~_mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                if (maskR == 0) {
                    while (true) {
                        if ((right - N) - left + 1 < 2 * N) {
                            goto end;
                        }

                        R = _mm256_loadu_si256((__m256i *) (array + right - N + 1));
                        const __m256i bytemask = _mm256_cmpgt_epi32(pivot, R);
                        if (_mm256_iszero(bytemask)) {
                            right -= N;
                        } else {
                            maskR = _mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                assert(left <= right);
                assert(maskL != 0);
                assert(maskR != 0);

                uint8_t mL;
                uint8_t mR;
                __m256i shuffleL;
                __m256i shuffleR;

                align_masks(maskL, maskR, mL, mR, shuffleL, shuffleR);
                swap_epi32(L, R, maskL, shuffleL, maskR, shuffleR);

                maskL = mL;
                maskR = mR;

                if (maskL == 0) {
                    _mm256_storeu_si256((__m256i *) (array + left), L);
                    left += N;
                }

                if (maskR == 0) {
                    _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
                    right -= N;
                }

            } // while

            end:

            assert(!(maskL != 0 && maskR != 0));

            if (maskL != 0) {
                _mm256_storeu_si256((__m256i *) (array + left), L);
            } else if (maskR != 0) {
                _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
            }

            if (left < right) {
                int less = 0;
                int greater = 0;
                const int all = right - left + 1;

                for (int i = left; i <= right; i++) {
                    less += int(array[i] < pv);
                    greater += int(array[i] > pv);
                }

                if (all == less) {
                    // all elements in range [left, right] less than pivot
                    scalar_partition_epi32(array, pv, origL, left);
                } else if (all == greater) {
                    // all elements in range [left, right] greater than pivot
                    scalar_partition_epi32(array, pv, left, origR);
                } else {
                    scalar_partition_epi32(array, pv, left, right);
                }
            }
        }

        //******************************************************************************
        //single precision floating point brojevi

        __m256 FORCE_INLINE bitmask_to_bytemask_ps(uint8_t bm) {

            const __m256 mask = _mm256_set1_ps(bm);
            const __m256 bits = _mm256_setr_ps(0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80);
            const __m256 tmp = _mm256_and_ps(mask, bits);

            return _mm256_cmp_ps(tmp, bits, _CMP_EQ_OQ);
        }


        void FORCE_INLINE
        align_masks_ps(uint8_t &a, uint8_t &b, uint8_t &rem_a, uint8_t &rem_b, __m256i &shuffle_a, __m256i &shuffle_b) {

            assert(a != 0);
            assert(b != 0);

            uint8_t tmpA = a;
            uint8_t tmpB = b;

            uint32_t __attribute__((__aligned__(32))) tmpshufa[8];
            uint32_t __attribute__((__aligned__(32))) tmpshufb[8];

            while (tmpA != 0 && tmpB != 0) {
                int idx_a = __builtin_ctz(tmpA);
                int idx_b = __builtin_ctz(tmpB);

                tmpA = tmpA & (tmpA - 1);
                tmpB = tmpB & (tmpB - 1);

                tmpshufb[idx_a] = idx_b;
                tmpshufa[idx_b] = idx_a;
            }

            a = a ^ tmpA;
            b = b ^ tmpB;

            assert(a != 0);
            assert(b != 0);
            assert(_mm_popcnt_u64(a) == _mm_popcnt_u64(b));

            rem_a = tmpA;
            rem_b = tmpB;

            shuffle_a = _mm256_load_si256((__m256i *) tmpshufa);
            shuffle_b = _mm256_load_si256((__m256i *) tmpshufb);
        }


        __m256 FORCE_INLINE merge_ps(const __m256 mask, const __m256 a, const __m256 b) {

            return _mm256_or_ps(
                    _mm256_and_ps(mask, a),
                    _mm256_andnot_ps(mask, b)
            );
        }


        void FORCE_INLINE swap_ps(
                __m256 &a, __m256 &b,
                uint8_t mask_a, const __m256i shuffle_a,
                uint8_t mask_b, const __m256i shuffle_b) {

            const __m256 to_swap_b = _mm256_permutevar8x32_ps(a, shuffle_a);
            const __m256 to_swap_a = _mm256_permutevar8x32_ps(b, shuffle_b);
            const __m256 ma = bitmask_to_bytemask_ps(mask_a);
            const __m256 mb = bitmask_to_bytemask_ps(mask_b);

            a = merge_ps(ma, to_swap_a, a);
            b = merge_ps(mb, to_swap_b, b);
        }


#define _mm256_iszero_ps(vec) (_mm256_testz_ps(vec, vec) != 0)

        void FORCE_INLINE partition_ps(float *array, float pv, int &left, int &right) {

            const int N = 8; // the number of items in a register (256/32)

            __m256 L;
            __m256 R;
            uint8_t maskL = 0;
            uint8_t maskR = 0;

            const __m256 pivot = _mm256_set1_ps(pv);

            int origL = left;
            int origR = right;

            while (true) {

                if (maskL == 0) {
                    while (true) {
                        if (right - (left + N) + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256*)(array + left)
                        L = _mm256_loadu_ps((const float *) (array + left));
                        const __m256 bytemask = _mm256_cmp_ps(pivot, L, _CMP_GT_OQ);

                        if (_mm256_testc_ps((__m256) bytemask, (__m256) _mm256_set1_ps(-1))) {
                            left += N;
                        } else {
                            maskL = ~_mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                if (maskR == 0) {
                    while (true) {
                        if ((right - N) - left + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256i *) (array + right - N + 1)
                        R = _mm256_loadu_ps((const float *) (array + right - N + 1));
                        const __m256 bytemask = _mm256_cmp_ps(pivot, R, _CMP_GT_OQ);
                        if (_mm256_iszero_ps(bytemask)) {
                            right -= N;
                        } else {
                            maskR = _mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                assert(left <= right);
                assert(maskL != 0);
                assert(maskR != 0);

                uint8_t mL;
                uint8_t mR;
                __m256i shuffleL;
                __m256i shuffleR;

                align_masks_ps(maskL, maskR, mL, mR, shuffleL, shuffleR);
                swap_ps(L, R, maskL, shuffleL, maskR, shuffleR);

                maskL = mL;
                maskR = mR;

                if (maskL == 0) {
                    _mm256_storeu_ps((float *) (array + left), L);
                    left += N;
                }

                if (maskR == 0) {
                    _mm256_storeu_ps((float *) (array + right - N + 1), R);
                    right -= N;
                }

            } // while

            end:

            assert(!(maskL != 0 && maskR != 0));

            if (maskL != 0) {
                _mm256_storeu_ps((float *) (array + left), L);
            } else if (maskR != 0) {
                _mm256_storeu_ps((float *) (array + right - N + 1), R);
            }

            if (left < right) {
                int less = 0;
                int greater = 0;
                const int all = right - left + 1;

                for (int i = left; i <= right; i++) {
                    less += int(array[i] < pv);
                    greater += int(array[i] > pv);
                }

                if (all == less) {
                    // all elements in range [left, right] less than pivot
                    scalar_partition_ps(array, pv, origL, left);
                } else if (all == greater) {
                    // all elements in range [left, right] greater than pivot
                    scalar_partition_ps(array, pv, left, origR);
                } else {
                    scalar_partition_ps(array, pv, left, right);
                }
            }
        }


        //****************************************************************************************
        //16 bitni integer brojevi

        __m256i FORCE_INLINE bitmask_to_bytemask_epi16(uint8_t bm) {

            const __m256i mask = _mm256_set1_epi16(bm);
            const __m256i bits = _mm256_setr_epi16(0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x100,0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x4000 );

            //const __m256i bits = _mm256_setr_epi16(0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08 );
            const __m256i tmp = _mm256_and_si256(mask, bits);

            return _mm256_cmpeq_epi16(tmp, bits);
        }


        void FORCE_INLINE align_masks_epi16(uint8_t &a, uint8_t &b, uint8_t &rem_a, uint8_t &rem_b, __m256i &shuffle_a, __m256i &shuffle_b) {

            assert(a != 0);
            assert(b != 0);

            uint8_t tmpA = a;
            uint8_t tmpB = b;

            uint32_t __attribute__((__aligned__(32))) tmpshufa[8];
            uint32_t __attribute__((__aligned__(32))) tmpshufb[8];

            while (tmpA != 0 && tmpB != 0) {
                int idx_a = __builtin_ctz(tmpA);
                int idx_b = __builtin_ctz(tmpB);

                tmpA = tmpA & (tmpA - 1);
                tmpB = tmpB & (tmpB - 1);

                tmpshufb[idx_a] = idx_b;
                tmpshufa[idx_b] = idx_a;
            }

            a = a ^ tmpA;
            b = b ^ tmpB;

            assert(a != 0);
            assert(b != 0);
            assert(_mm_popcnt_u64(a) == _mm_popcnt_u64(b));

            rem_a = tmpA;
            rem_b = tmpB;

            shuffle_a = _mm256_load_si256((__m256i *) tmpshufa);
            shuffle_b = _mm256_load_si256((__m256i *) tmpshufb);
        }


        __m256i FORCE_INLINE merge_epi16(const __m256i mask, const __m256i a, const __m256i b) {

            return _mm256_or_si256(
                    _mm256_andnot_si256(mask, a),
                    _mm256_andnot_si256(mask, b)
            );
        }


        void FORCE_INLINE swap_epi16(
                __m256i &a, __m256i &b,
                uint8_t mask_a, const __m256i shuffle_a,
                uint8_t mask_b, const __m256i shuffle_b) {

            const __m256i to_swap_b = _mm256_permutevar8x32_epi32(a, shuffle_a);
            const __m256i to_swap_a = _mm256_permutevar8x32_epi32(b, shuffle_b);
            const __m256i ma = bitmask_to_bytemask_epi16(mask_a);
            const __m256i mb = bitmask_to_bytemask_epi16(mask_b);

            a = merge_epi16(ma, to_swap_a, a);
            b = merge_epi16(mb, to_swap_b, b);
        }


#define _mm256_iszero(vec) (_mm256_testz_si256(vec, vec) != 0)

        void FORCE_INLINE partition_epi16(uint16_t *array, uint16_t pv, int &left, int &right) {

            const int N = 16; // the number of items in a register (256/16)

            __m256i L;
            __m256i R;
            uint8_t maskL = 0;
            uint8_t maskR = 0;

            const __m256i pivot = _mm256_set1_epi16(pv);

            int origL = left;
            int origR = right;

            while (true) {

                if (maskL == 0) {
                    while (true) {
                        if (right - (left + N) + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256*)(array + left)
                        L = _mm256_loadu_si256((__m256i *) (array + left));
                        const __m256i bytemask = _mm256_cmpgt_epi16(pivot, L);

                        if (_mm256_testc_ps((__m256) bytemask, (__m256) _mm256_set1_epi16(-1))) {
                            left += N;
                        } else {
                            maskL = ~_mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                if (maskR == 0) {
                    while (true) {
                        if ((right - N) - left + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256i *) (array + right - N + 1)
                        R = _mm256_loadu_si256((__m256i *) (array + right - N + 1));
                        const __m256i bytemask = _mm256_cmpgt_epi16(pivot, R);
                        if (_mm256_iszero(bytemask)) {
                            right -= N;
                        } else {
                            maskR = _mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                assert(left <= right);
                assert(maskL != 0);
                assert(maskR != 0);

                uint8_t mL;
                uint8_t mR;
                __m256i shuffleL;
                __m256i shuffleR;

                align_masks_epi16(maskL, maskR, mL, mR, shuffleL, shuffleR);
                swap_epi16(L, R, maskL, shuffleL, maskR, shuffleR);

                maskL = mL;
                maskR = mR;

                if (maskL == 0) {
                    _mm256_storeu_si256((__m256i *) (array + left), L);
                    left += N;
                }

                if (maskR == 0) {
                    _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
                    right -= N;
                }

            } // while

            end:

            assert(!(maskL != 0 && maskR != 0));

            if (maskL != 0) {
                _mm256_storeu_si256((__m256i *) (array + left), L);
            } else if (maskR != 0) {
                _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
            }

            if (left < right) {
                int less = 0;
                int greater = 0;
                const int all = right - left + 1;

                for (int i = left; i <= right; i++) {
                    less += int(array[i] < pv);
                    greater += int(array[i] > pv);
                }

                if (all == less) {
                    // all elements in range [left, right] less than pivot
                    scalar_partition_16(array, pv, origL, left);
                } else if (all == greater) {
                    // all elements in range [left, right] greater than pivot
                    scalar_partition_16(array, pv, left, origR);
                } else {
                    scalar_partition_16(array, pv, left, right);
                }
            }
        }
        //***************************************************************************
        //za 64-bitne integer brojeve

        __m256i FORCE_INLINE bitmask_to_bytemask_epi64(uint8_t bm) {

            const __m256i mask = _mm256_set1_epi64x(bm);
            const __m256i bits = _mm256_setr_epi64x(0x01, 0x02, 0x04, 0x08 );
            const __m256i tmp = _mm256_and_si256(mask, bits);

            return _mm256_cmpeq_epi64(tmp, bits);
        }


        void FORCE_INLINE align_masks_epi64(uint8_t &a, uint8_t &b, uint8_t &rem_a, uint8_t &rem_b, __m256i &shuffle_a, __m256i &shuffle_b) {

            assert(a != 0);
            assert(b != 0);

            uint8_t tmpA = a;
            uint8_t tmpB = b;

            uint32_t __attribute__((__aligned__(32))) tmpshufa[8];
            uint32_t __attribute__((__aligned__(32))) tmpshufb[8];

            while (tmpA != 0 && tmpB != 0) {
                int idx_a = __builtin_ctz(tmpA);
                int idx_b = __builtin_ctz(tmpB);

                tmpA = tmpA & (tmpA - 1);
                tmpB = tmpB & (tmpB - 1);

                tmpshufb[idx_a] = idx_b;
                tmpshufa[idx_b] = idx_a;
            }

            a = a ^ tmpA;
            b = b ^ tmpB;

            assert(a != 0);
            assert(b != 0);
            assert(_mm_popcnt_u64(a) == _mm_popcnt_u64(b));

            rem_a = tmpA;
            rem_b = tmpB;

            shuffle_a = _mm256_load_si256((__m256i *) tmpshufa);
            shuffle_b = _mm256_load_si256((__m256i *) tmpshufb);
        }


        __m256i FORCE_INLINE merge_epi64(const __m256i mask, const __m256i a, const __m256i b) {

            return _mm256_or_si256(
                    _mm256_andnot_si256(mask, a),
                    _mm256_andnot_si256(mask, b)
            );
        }


        void FORCE_INLINE swap_epi64(
                __m256i &a, __m256i &b,
                uint8_t mask_a, const __m256i shuffle_a,
                uint8_t mask_b, const __m256i shuffle_b) {

            const __m256i to_swap_b = _mm256_permutevar8x32_epi32(a, shuffle_a);
            const __m256i to_swap_a = _mm256_permutevar8x32_epi32(b, shuffle_b);
            const __m256i ma = bitmask_to_bytemask_epi64(mask_a);
            const __m256i mb = bitmask_to_bytemask_epi64(mask_b);

            a = merge_epi64(ma, to_swap_a, a);
            b = merge_epi64(mb, to_swap_b, b);
        }


#define _mm256_iszero(vec) (_mm256_testz_si256(vec, vec) != 0)

        void FORCE_INLINE partition_epi64(uint64_t *array, uint64_t pv, int &left, int &right) {

            const int N = 4; // the number of items in a register (256/64)

            __m256i L;
            __m256i R;
            uint8_t maskL = 0;
            uint8_t maskR = 0;

            const __m256i pivot = _mm256_set1_epi64x(pv);

            int origL = left;
            int origR = right;

            while (true) {

                if (maskL == 0) {
                    while (true) {
                        if (right - (left + N) + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256*)(array + left)
                        L = _mm256_loadu_si256((__m256i *) (array + left));
                        const __m256i bytemask = _mm256_cmpgt_epi64(pivot, L);

                        if (_mm256_testc_ps((__m256) bytemask, (__m256) _mm256_set1_epi64x(-1))) {
                            left += N;
                        } else {
                            maskL = ~_mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                if (maskR == 0) {
                    while (true) {
                        if ((right - N) - left + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256i *) (array + right - N + 1)
                        R = _mm256_loadu_si256((__m256i *) (array + right - N + 1));
                        const __m256i bytemask = _mm256_cmpgt_epi64(pivot, R);
                        if (_mm256_iszero(bytemask)) {
                            right -= N;
                        } else {
                            maskR = _mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                assert(left <= right);
                assert(maskL != 0);
                assert(maskR != 0);

                uint8_t mL;
                uint8_t mR;
                __m256i shuffleL;
                __m256i shuffleR;

                align_masks_epi64(maskL, maskR, mL, mR, shuffleL, shuffleR);
                swap_epi64(L, R, maskL, shuffleL, maskR, shuffleR);

                maskL = mL;
                maskR = mR;

                if (maskL == 0) {
                    _mm256_storeu_si256((__m256i *) (array + left), L);
                    left += N;
                }

                if (maskR == 0) {
                    _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
                    right -= N;
                }

            } // while

            end:

            assert(!(maskL != 0 && maskR != 0));

            if (maskL != 0) {
                _mm256_storeu_si256((__m256i *) (array + left), L);
            } else if (maskR != 0) {
                _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
            }

            if (left < right) {
                int less = 0;
                int greater = 0;
                const int all = right - left + 1;

                for (int i = left; i <= right; i++) {
                    less += int(array[i] < pv);
                    greater += int(array[i] > pv);
                }

                if (all == less) {
                    // all elements in range [left, right] less than pivot
                    scalar_partition_64(array, pv, origL, left);
                } else if (all == greater) {
                    // all elements in range [left, right] greater than pivot
                    scalar_partition_64(array, pv, left, origR);
                } else {
                    scalar_partition_64(array, pv, left, right);
                }
            }
        }

        //************************************************************************************************
        //za 8bitne integer brojeve

        __m256i FORCE_INLINE bitmask_to_bytemask_epi8(uint8_t bm) {

            const __m256i mask = _mm256_set1_epi8(bm);
            //const __m256i bits = _mm256_setr_epi8(0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x100,0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000,
            //                                      0x10000, 0x20000, 0x40000, 0x80000, 0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, 0x4000000, 0x8000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000);
            const __m256i bits = _mm256_setr_epi8(0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08,
                                                  0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08, 0x01, 0x02, 0x04, 0x08);

            const __m256i tmp = _mm256_and_si256(mask, bits);

            return _mm256_cmpeq_epi8(tmp, bits);
        }


        void FORCE_INLINE align_masks_epi8(uint8_t &a, uint8_t &b, uint8_t &rem_a, uint8_t &rem_b, __m256i &shuffle_a, __m256i &shuffle_b) {

            assert(a != 0);
            assert(b != 0);

            uint8_t tmpA = a;
            uint8_t tmpB = b;

            uint32_t __attribute__((__aligned__(32))) tmpshufa[8];
            uint32_t __attribute__((__aligned__(32))) tmpshufb[8];

            while (tmpA != 0 && tmpB != 0) {
                int idx_a = __builtin_ctz(tmpA);
                int idx_b = __builtin_ctz(tmpB);

                tmpA = tmpA & (tmpA - 1);
                tmpB = tmpB & (tmpB - 1);

                tmpshufb[idx_a] = idx_b;
                tmpshufa[idx_b] = idx_a;
            }

            a = a ^ tmpA;
            b = b ^ tmpB;

            assert(a != 0);
            assert(b != 0);
            assert(_mm_popcnt_u64(a) == _mm_popcnt_u64(b));

            rem_a = tmpA;
            rem_b = tmpB;

            shuffle_a = _mm256_load_si256((__m256i *) tmpshufa);
            shuffle_b = _mm256_load_si256((__m256i *) tmpshufb);
        }


        __m256i FORCE_INLINE merge_epi8(const __m256i mask, const __m256i a, const __m256i b) {

            return _mm256_or_si256(
                    _mm256_andnot_si256(mask, a),
                    _mm256_andnot_si256(mask, b)
            );
        }


        void FORCE_INLINE swap_epi8(
                __m256i &a, __m256i &b,
                uint8_t mask_a, const __m256i shuffle_a,
                uint8_t mask_b, const __m256i shuffle_b) {

            const __m256i to_swap_b = _mm256_permutevar8x32_epi32(a, shuffle_a);
            const __m256i to_swap_a = _mm256_permutevar8x32_epi32(b, shuffle_b);
            const __m256i ma = bitmask_to_bytemask_epi8(mask_a);
            const __m256i mb = bitmask_to_bytemask_epi8(mask_b);

            a = merge_epi8(ma, to_swap_a, a);
            b = merge_epi8(mb, to_swap_b, b);
        }


#define _mm256_iszero(vec) (_mm256_testz_si256(vec, vec) != 0)

        void FORCE_INLINE partition_epi8(uint8_t *array, uint8_t pv, int &left, int &right) {

            const int N = 32; // the number of items in a register (256/8)

            __m256i L;
            __m256i R;
            uint8_t maskL = 0;
            uint8_t maskR = 0;

            const __m256i pivot = _mm256_set1_epi8(pv);

            int origL = left;
            int origR = right;

            while (true) {

                if (maskL == 0) {
                    while (true) {
                        if (right - (left + N) + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256*)(array + left)
                        L = _mm256_loadu_si256((__m256i *) (array + left));
                        const __m256i bytemask = _mm256_cmpgt_epi8(pivot, L);

                        if (_mm256_testc_ps((__m256) bytemask, (__m256) _mm256_set1_epi8(-1))) {
                            left += N;
                        } else {
                            maskL = ~_mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                if (maskR == 0) {
                    while (true) {
                        if ((right - N) - left + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256i *) (array + right - N + 1)
                        R = _mm256_loadu_si256((__m256i *) (array + right - N + 1));
                        const __m256i bytemask = _mm256_cmpgt_epi8(pivot, R);
                        if (_mm256_iszero(bytemask)) {
                            right -= N;
                        } else {
                            maskR = _mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                assert(left <= right);
                assert(maskL != 0);
                assert(maskR != 0);

                uint8_t mL;
                uint8_t mR;
                __m256i shuffleL;
                __m256i shuffleR;

                align_masks_epi8(maskL, maskR, mL, mR, shuffleL, shuffleR);
                swap_epi8(L, R, maskL, shuffleL, maskR, shuffleR);

                maskL = mL;
                maskR = mR;

                if (maskL == 0) {
                    _mm256_storeu_si256((__m256i *) (array + left), L);
                    left += N;
                }

                if (maskR == 0) {
                    _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
                    right -= N;
                }

            } // while

            end:

            assert(!(maskL != 0 && maskR != 0));

            if (maskL != 0) {
                _mm256_storeu_si256((__m256i *) (array + left), L);
            } else if (maskR != 0) {
                _mm256_storeu_si256((__m256i *) (array + right - N + 1), R);
            }

            if (left < right) {
                int less = 0;
                int greater = 0;
                const int all = right - left + 1;

                for (int i = left; i <= right; i++) {
                    less += int(array[i] < pv);
                    greater += int(array[i] > pv);
                }

                if (all == less) {
                    // all elements in range [left, right] less than pivot
                    scalar_partition_8(array, pv, origL, left);
                } else if (all == greater) {
                    // all elements in range [left, right] greater than pivot
                    scalar_partition_8(array, pv, left, origR);
                } else {
                    scalar_partition_8(array, pv, left, right);
                }
            }
        }

        //*********************************************************************************************************
        //za double
        __m256d FORCE_INLINE bitmask_to_bytemask_pd(uint8_t bm) {

            const __m256d mask = _mm256_set1_pd(bm);
            const __m256d bits = _mm256_setr_pd(0x01, 0x02, 0x04, 0x08);
            const __m256d tmp = _mm256_and_pd(mask, bits);

            return _mm256_cmp_pd(tmp, bits, _CMP_EQ_OQ);
        }


        void FORCE_INLINE align_masks_pd(uint8_t &a, uint8_t &b, uint8_t &rem_a, uint8_t &rem_b, __m256i &shuffle_a, __m256i &shuffle_b) {

            assert(a != 0);
            assert(b != 0);

            uint8_t tmpA = a;
            uint8_t tmpB = b;

            uint32_t __attribute__((__aligned__(32))) tmpshufa[8];
            uint32_t __attribute__((__aligned__(32))) tmpshufb[8];

            while (tmpA != 0 && tmpB != 0) {
                int idx_a = __builtin_ctz(tmpA);
                int idx_b = __builtin_ctz(tmpB);

                tmpA = tmpA & (tmpA - 1);
                tmpB = tmpB & (tmpB - 1);

                tmpshufb[idx_a] = idx_b;
                tmpshufa[idx_b] = idx_a;
            }

            a = a ^ tmpA;
            b = b ^ tmpB;

            assert(a != 0);
            assert(b != 0);
            assert(_mm_popcnt_u64(a) == _mm_popcnt_u64(b));

            rem_a = tmpA;
            rem_b = tmpB;

            shuffle_a = _mm256_load_si256((__m256i *) tmpshufa);
            shuffle_b = _mm256_load_si256((__m256i *) tmpshufb);
        }


        __m256d FORCE_INLINE merge_pd(const __m256d mask, const __m256d a, const __m256d b) {

            return _mm256_or_pd(
                    _mm256_and_pd(mask, a),
                    _mm256_andnot_pd(mask, b)
            );
        }


        void FORCE_INLINE swap_pd(
                __m256d &a, __m256d &b,
                uint8_t mask_a, const __m256i shuffle_a,
                uint8_t mask_b, const __m256i shuffle_b) {

            const __m256d to_swap_b = _mm256_permutevar_pd(a, shuffle_a);
            const __m256d to_swap_a = _mm256_permutevar_pd(b, shuffle_b);
            const __m256d ma = bitmask_to_bytemask_pd(mask_a);
            const __m256d mb = bitmask_to_bytemask_pd(mask_b);

            a = merge_pd(ma, to_swap_a, a);
            b = merge_pd(mb, to_swap_b, b);
        }


#define _mm256_iszero_pd(vec) (_mm256_testz_pd(vec, vec) != 0)

        void FORCE_INLINE partition_pd(double *array, double pv, int &left, int &right) {

            const int N = 4; // the number of items in a register (256/64)

            __m256d L;
            __m256d R;
            uint8_t maskL = 0;
            uint8_t maskR = 0;

            const __m256d pivot = _mm256_set1_pd(pv);

            int origL = left;
            int origR = right;

            while (true) {

                if (maskL == 0) {
                    while (true) {
                        if (right - (left + N) + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256*)(array + left)
                        L = _mm256_loadu_pd((const double *) (array + left));
                        const __m256d bytemask = _mm256_cmp_pd(pivot, L, _CMP_GT_OQ);

                        if (_mm256_testc_ps((__m256) bytemask, (__m256) _mm256_set1_pd(-1))) {
                            left += N;
                        } else {
                            maskL = ~_mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                if (maskR == 0) {
                    while (true) {
                        if ((right - N) - left + 1 < 2 * N) {
                            goto end;
                        }
                        // (__m256i *) (array + right - N + 1)
                        R = _mm256_loadu_pd((const double *) (array + right - N + 1));
                        const __m256d bytemask = _mm256_cmp_pd(pivot, R, _CMP_GT_OQ);
                        if (_mm256_iszero_pd(bytemask)) {
                            right -= N;
                        } else {
                            maskR = _mm256_movemask_ps((__m256) bytemask);
                            break;
                        }
                    }

                }

                assert(left <= right);
                assert(maskL != 0);
                assert(maskR != 0);

                uint8_t mL;
                uint8_t mR;
                __m256i shuffleL;
                __m256i shuffleR;

                align_masks_pd(maskL, maskR, mL, mR, shuffleL, shuffleR);
                swap_pd(L, R, maskL, shuffleL, maskR, shuffleR);

                maskL = mL;
                maskR = mR;

                if (maskL == 0) {
                    _mm256_storeu_pd((double *) (array + left), L);
                    left += N;
                }

                if (maskR == 0) {
                    _mm256_storeu_pd((double *) (array + right - N + 1), R);
                    right -= N;
                }

            } // while

            end:

            assert(!(maskL != 0 && maskR != 0));

            if (maskL != 0) {
                _mm256_storeu_pd((double *) (array + left), L);
            } else if (maskR != 0) {
                _mm256_storeu_pd((double *) (array + right - N + 1), R);
            }

            if (left < right) {
                int less = 0;
                int greater = 0;
                const int all = right - left + 1;

                for (int i = left; i <= right; i++) {
                    less += int(array[i] < pv);
                    greater += int(array[i] > pv);
                }

                if (all == less) {
                    // all elements in range [left, right] less than pivot
                    scalar_partition_pd(array, pv, origL, left);
                } else if (all == greater) {
                    // all elements in range [left, right] greater than pivot
                    scalar_partition_pd(array, pv, left, origR);
                } else {
                    scalar_partition_pd(array, pv, left, right);
                }
            }
        }
    }// namespace avx2

} // namespace qs

