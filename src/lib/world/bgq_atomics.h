/*
 * Written by Jeff Hammond, July 2012
 * Copyright Argonne National Laboratory
 *
 * This implementation of BGQ atomics is based upon 
 * hwi/include/bqc/A2_inlines.h but uses signed integers
 * instead of unsigned integer and/or long types.
 */

/* for ppc_msync only; should just implemente here and omit header inclusion */
#include <hwi/include/bqc/A2_inlines.h>

__INLINE__ int32_t LoadReservedSigned32( volatile int32_t *pVar )
{
   register int32_t Val;
   asm volatile ("lwarx   %[rc],0,%[pVar];"
                 : [rc] "=&b" (Val)
                 : [pVar] "b" (pVar));
   return(Val);
}

__INLINE__ int StoreConditionalSigned32( volatile int32_t *pVar, int32_t Val )
{
   register int rc = 1; // assume success
   asm volatile ("  stwcx.  %2,0,%1;"
                 "  beq     1f;"       // conditional store succeeded
                 "  li      %0,0;"
                 "1:;"
                 : "=b" (rc)
                 : "b"  (pVar),
                   "b"  (Val),
                   "0"  (rc)
                 : "cc", "memory" );
   return(rc);
}

__INLINE__ int32_t CompareAndSwapSigned32( volatile int32_t *var,
                                           int32_t  Compare,
                                           int32_t  NewValue )
{
    int32_t OldValue = *var;

    do {
       int32_t TmpValue = LoadReservedSigned32( var );
       if ( Compare != TmpValue  )
          {
              return(OldValue);
          }
       }
       while( !StoreConditionalSigned32( var, NewValue ) );

    return(OldValue);
}

__INLINE__ int32_t FetchAndAddSigned32( volatile int32_t *pVar, int32_t value )
{
    register int32_t old_val, tmp_val;

    ppc_msync();

    do
    {
        old_val = LoadReservedSigned32( pVar );
        tmp_val = old_val + value;
    }
    while ( !StoreConditionalSigned32( pVar, tmp_val ) );

    return( old_val );
}
