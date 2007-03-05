# Copyright (c) 2006, R.J. Harrison, UT/ORNL
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the names of Oak Ridge National Laboratory and the
#       University of Tennessee, nor the names of its contributors may
#       be used to endorse or promote products derived from this
#       software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

        
        
        
       	.file	"ppcaix32.s"
	.toc
	.csect .text[PR]
	.align 2
	.globl fiber_make_context
	.globl .fiber_make_context
	.csect fiber_make_context[DS]
fiber_make_context:
	.long .fiber_make_context, TOC[tc0], 0
	.csect .text[PR]
.fiber_make_context:
	mflr 0

	# Arguments
	# r3 is new stack pointer
	# r4 is function pointer (toc entry)

	# Return value
	# updated stack pointer
	
        # Saved stack is (with all data 8 bytes for easy maint. of 32/64 versions)
	# cr
	# r0 (return address if any)
	# address to resume execution
	# toc (of the to be called function)
	# 18 saved gprs 
	# 18 saved fprs
	
	mfcr 5
	stw  5, -8(3)   # cr
	stw  0, -16(3)  # return address
	lwz  5, 0(4)
	stw  5, -24(3)  # address of routine to invoke
	lwz  5, 4(4)
	stw  5, -32(3)  # new routine toc

        addi 3,3,-32

	stw  14,-144(3)	
	stw  15,-136(3)	
	stw  16,-128(3)	
	stw  17,-120(3)	
	stw  18,-112(3)	
	stw  19,-104(3)	
	stw  20,-96(3) 	
	stw  21,-88(3) 	
	stw  22,-80(3) 	
	stw  23,-72(3) 	
	stw  24,-64(3) 	
	stw  25,-56(3) 	
	stw  26,-48(3) 	
	stw  27,-40(3) 	
	stw  28,-32(3) 	
	stw  29,-24(3) 	
	stw  30,-16(3) 	
	stw  31,-8(3)

	addi 3,3,-144
			
	stfd 14,-144(3) 
	stfd 15,-136(3) 	
	stfd 16,-128(3) 	
	stfd 17,-120(3) 	
	stfd 18,-112(3) 	
	stfd 19,-104(3) 	
	stfd 20,-96(3)  	
	stfd 21,-88(3)  	
	stfd 22,-80(3)  	
	stfd 23,-72(3)  	
	stfd 24,-64(3)  	
	stfd 25,-56(3)  		
	stfd 26,-48(3)  	
	stfd 27,-40(3)  	
	stfd 28,-32(3)  	
	stfd 29,-24(3)  	
	stfd 30,-16(3)  	
	stfd 31,-8(3)   
		
	addi 3, 3, -144

	mtlr 0
	blr
	nop

	.toc
	.csect .text[PR]
	.align 2
	.globl fiber_set_context
	.globl .fiber_set_context
	.csect fiber_set_context[DS]
fiber_set_context:
	.long .fiber_set_context, TOC[tc0], 0
	.csect .text[PR]
.fiber_set_context:
	mflr 0

	# Arguments
	# r3 is new stack pointer

	# Return value
	# Never returns
	
        # Saved stack is (with all data 8 bytes for easy maint. of 32/64 versions)
	# cr
	# r0 (return address if any)
	# address to resume execution
	# toc (of the to be called function)
	# 18 saved gprs 
	# 18 saved fprs
	
        addi 3, 3, 144

	lfd 14,-144(3) 
	lfd 15,-136(3) 	
	lfd 16,-128(3) 	
	lfd 17,-120(3) 	
	lfd 18,-112(3) 	
	lfd 19,-104(3) 	
	lfd 20,-96(3)  	
	lfd 21,-88(3)  	
	lfd 22,-80(3)  	
	lfd 23,-72(3)  	
	lfd 24,-64(3)  	
	lfd 25,-56(3)  		
	lfd 26,-48(3)  	
	lfd 27,-40(3)  	
	lfd 28,-32(3)  	
	lfd 29,-24(3)  	
	lfd 30,-16(3)  	
	lfd 31,-8(3)   

	addi 3, 3, 144
		
	lwz  14,-144(3)	
	lwz  15,-136(3)	
	lwz  16,-128(3)	
	lwz  17,-120(3)	
	lwz  18,-112(3)	
	lwz  19,-104(3)	
	lwz  20,-96(3) 	
	lwz  21,-88(3) 	
	lwz  22,-80(3) 	
	lwz  23,-72(3) 	
	lwz  24,-64(3) 	
	lwz  25,-56(3) 	
	lwz  26,-48(3) 	
	lwz  27,-40(3) 	
	lwz  28,-32(3) 	
	lwz  29,-24(3) 	
	lwz  30,-16(3) 	
	lwz  31,-8(3)

	addi 3, 3, 32	

	lwz 2, -32(3)    # restore old toc
	lwz 4, -24(3)    # address to resume execution
	lwz 0, -16(3)    # restore old return address
	lwz 5, -8(3)     # restore old cr
	mtcr 5
	
	mr 1,3           # restore old stack

	mtctr 4
	bctr

	
	.file	"ppcaix32.s"
	.toc
	.csect .text[PR]
	.align 2
	.globl fiber_swap_context
	.globl .fiber_swap_context
	.csect fiber_swap_context[DS]
fiber_swap_context:
	.long .fiber_swap_context, TOC[tc0], 0
	.csect .text[PR]
.fiber_swap_context:
	mflr 0

	# Arguments
	# r3 points to location to store old stack pointer with saved state
	# r4 is new stack pointer

	# Return value
	# None (returns upon fiber resume)
	
        # Saved stack is (with all data 8 bytes for easy maint. of 32/64 versions)
	# cr
	# new stack ptr (everything saved will be popped)
	# address to resume execution
	# toc (of the to be called function)
	# 18 saved gprs 
	# 18 saved fprs
	
	mfcr 5
	stw  5, -8(1)   # cr
	stw  0, -16(1)  # old return address
	stw  0, -24(1)  # address to which control will return on resume
	stw  2, -32(1)  # toc

        addi 1,1,-32

	stw  14,-144(1)	
	stw  15,-136(1)	
	stw  16,-128(1)	
	stw  17,-120(1)	
	stw  18,-112(1)	
	stw  19,-104(1)	
	stw  20,-96(1) 	
	stw  21,-88(1) 	
	stw  22,-80(1) 	
	stw  23,-72(1) 	
	stw  24,-64(1) 	
	stw  25,-56(1) 	
	stw  26,-48(1) 	
	stw  27,-40(1) 	
	stw  28,-32(1) 	
	stw  29,-24(1) 	
	stw  30,-16(1) 	
	stw  31,-8(1)

	addi 1,1,-144
			
	stfd 14,-144(1) 
	stfd 15,-136(1) 	
	stfd 16,-128(1) 	
	stfd 17,-120(1) 	
	stfd 18,-112(1) 	
	stfd 19,-104(1) 	
	stfd 20,-96(1)  	
	stfd 21,-88(1)  	
	stfd 22,-80(1)  	
	stfd 23,-72(1)  	
	stfd 24,-64(1)  	
	stfd 25,-56(1)  		
	stfd 26,-48(1)  	
	stfd 27,-40(1)  	
	stfd 28,-32(1)  	
	stfd 29,-24(1)  	
	stfd 30,-16(1)  	
	stfd 31,-8(1)    
		
	addi 1, 1, -144

	stw 1, 0(3)  # Save the old stack pointer with saved state

	mr 3, 4    # Move new stack ptr (4) into first argument register (3) for set_context
	b .fiber_set_context 

	
	
