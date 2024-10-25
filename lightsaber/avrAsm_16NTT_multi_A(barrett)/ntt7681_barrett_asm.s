
/*!
 * @details     void q1_barrett_ntt_asm (int16_t a[256]);
 * @brief       forword NTT 7681 
 * @param[in]   signed 16-bit x 256 polynomial data : a
 * @param[out]  return NTT(a)
 * @note        ..
 * @pre         ..
 * @author Youngbeom Kim
 */
            EXTERN  q2_streamlined_CT_table
            NAME    q2_barrett_ntt_asm      ;   module name
            PUBLIC  q2_barrett_ntt_asm      ;   make the main label vissible outside this module
            RSEG    CODE                     ;   place program in 'CODE' segment
    

mc_muls16x16_32    MACRO    m_a0, m_a1, m_b0, m_b1, m_c0, m_c1, m_c2, m_c3, m_zero

	muls    m_a1,   m_b1        
	movw	m_c2,   r0
	mul     m_a0,   m_b0        
	movw	m_c0,   r0
	mulsu	m_a1,   m_b0        
	sbc     m_c3,   m_zero
	add     m_c1,   r0
	adc     m_c2,   r1
	adc     m_c3,   m_zero
	mulsu	m_b1,   m_a0        
	sbc     m_c3,   m_zero
	add     m_c1,   r0
	adc     m_c2,   r1
	adc     m_c3,   m_zero

    ENDM
    
mc_muls16x16_lo16    MACRO    m_a0, m_a1, m_b0, m_b1, m_c0, m_c1, m_zero

	mul     m_a0,   m_b0
    movw    m_c0,   r0
    mulsu   m_a1,   m_b0
    add     m_c1,   r0
    mulsu   m_b1,   m_a0
    add     m_c1,   r0
    
    ENDM    
    
mc_ct_hi        MACRO    m_a0, m_a1, m_b0, m_b1

	ld      m_a0,   X+
    ld      m_a1,   X
    ldd     m_b0,   Y+0
    ldd     m_b1,   Y+1

    ENDM
    
mc_barrett_mul  MACRO    m_b0, m_b1, m_z0, m_z1, m_bp0, m_bp1, m_lo0, m_lo1, m_brt0, m_brt1, m_brt2, m_brt3, m_c_tmp, m_cst_q, m_zero

    mc_muls16x16_lo16   m_b0, m_b1, m_z0,  m_z1,  m_lo0,  m_lo1,  m_zero
    mc_muls16x16_32     m_b0, m_b1, m_bp0, m_bp1, m_brt0, m_brt1, m_brt2, m_brt3, m_zero
    ldi     c_tmp,      0x80
    add     m_brt1,     m_c_tmp
    adc     m_brt2,     m_zero
    adc     m_brt3,     m_zero    
    movw    m_b0,       m_brt2
    ldi     m_cst_q,    0x1E ; 7681
    mul     m_b0,       m_cst_q
    add     m_b1,       r0    
    sub     m_lo0,      m_b0
    sbc     m_lo1,      m_b1    
    
    ENDM

mc_ct_lo        MACRO    m_a0, m_a1, m_b0, m_b1, m_tmp0, m_tmp1

	movw    m_tmp0, m_a0
    add     m_a0,   m_b0
    adc     m_a1,   m_b1
    sub     m_tmp0, m_b0
    sbc     m_tmp1, m_b1    
    
    st      X,      m_a1
    st      -X,     m_a0
    st      Y+,     m_tmp0
    st      Y+,     m_tmp1    
    adiw    r26,    2

    ENDM    
                      
q2_barrett_ntt_asm:
    ;prolog
    push    r0
    push    r1
    push    r2
    push    r3
    push    r4
    push    r5
    push    r6
    push    r7
    push    r8
    push    r9
    push    r10
    push    r11
    push    r12
    push    r13
    push    r14
    push    r15
    push    r16
    push    r17
    push    r18
    push    r19
    push    r20
    push    r21
    push    r22
    push    r23
    push    r24
    push    r25
    push    r26
    push    r27
    push    r28
    push    r29
    push    r30
    push    r31

#define	zero        r2
#define tmp0        r3
#define tmp1        r4
#define	offset 	    r5	
#define	a0     	    r6	
#define	a1     	    r7	
#define	adr0   	    r8	
#define	adr1   	    r9	
#define	brt0        r10	
#define	brt1        r11	
#define	brt2        r12	
#define	brt3        r13	
#define	lo0    	    r14	
#define	lo1    	    r15	

#define	cst_q  	    r16	
#define	c_tmp	    r17	
#define	a_t0  	    r16	
#define	a_t1	    r17	

#define	b0     	    r18	
#define	b1     	    r19	
#define	bp0    	    r20	
#define	bp1    	    r21	
#define	z0     	    r22	
#define	z1     	    r23	
#define	i_cnt  	    r24	
#define	o_cnt  	    r25

#define	Xlo  	    r26
#define	Xhi  	    r27
#define	Ylo  	    r28
#define	Yhi  	    r29
#define	Zlo  	    r30
#define	Zhi  	    r31


initialize:
    clr     zero
    ;   X = (int16_t*)&a, Y= (int16_t*)&a + 128
    movw    Xlo,    r16
    movw    Ylo,    r16    
    movw    adr0,   r16
    
    ;   Z = zetas || zetashi table
    ldi     r30,    low (q2_streamlined_CT_table)
    ldi     r31,    high(q2_streamlined_CT_table)    
    
layer0: 
    inc     r29 ; layer 0 addr

    ldi     i_cnt,  128   
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+ 
        
in0:    
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    dec     i_cnt
    brne    in0
    
    
layer1: 
    ldi     i_cnt,  64      ;
    ldi     o_cnt,  2       ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  128     ;
    mov     offset, c_tmp
    
    add     r28,    c_tmp
    adc     r29,    zero
    
out1:        
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+     
in1:
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1   
    dec     i_cnt
    brne    in1             ;
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    ldi     i_cnt,  64      ;
    
    dec     o_cnt
    brne    out1            ;


layer2: 
    ldi     i_cnt,  32      ;
    ldi     o_cnt,  4       ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  64     ;
    mov     offset, c_tmp
    
    add     r28,    c_tmp
    adc     r29,    zero
    
out2:        
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+      
in2:
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    
    dec     i_cnt
    brne    in2             ;
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    ldi     i_cnt,  32      ;
    
    dec     o_cnt
    brne    out2            ;


layer3:  
    ldi     i_cnt,  16      ;
    ldi     o_cnt,  8       ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  32      ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
    
out3:        
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+     
in3:
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    
    dec     i_cnt
    brne    in3             ;
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    ldi     i_cnt,  16      ;
    
    dec     o_cnt
    brne    out3            ;


layer4: 
    ldi     i_cnt,  8       ;
    ldi     o_cnt,  16      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  16      ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
    
out4:        
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+   
in4:
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1   
    
    dec     i_cnt
    brne    in4             ;
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    ldi     i_cnt,  8      ;
    
    dec     o_cnt
    brne    out4            ;

layer5:  
    ldi     i_cnt,  4       ;
    ldi     o_cnt,  32      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  8       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
    
out5:        
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+  
in5:
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1  
    
    dec     i_cnt
    brne    in5             ;
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    ldi     i_cnt,  4      ;
    
    dec     o_cnt
    brne    out5            ;

layer6: 
    ldi     i_cnt,  2       ;
    ldi     o_cnt,  64      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  4       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
    
out6:        
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+   
in6:
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    
    dec     i_cnt
    brne    in6             ;
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    ldi     i_cnt,  2      ;
    
    dec     o_cnt
    breq    layer7          ;
    rjmp    out6
    
layer7: // Layer 6
    ldi     o_cnt,  128      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  2       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
    
out7:        
    lpm     z0,     Z+
    lpm     z1,     Z+  
    lpm     bp0,    Z+
    lpm     bp1,    Z+ 
    
    mc_ct_hi        a0, a1, b0, b1
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1  
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    dec     o_cnt
    brne    out7    
    
end_ntt:

#undef	zero        
#undef  tmp0        
#undef  tmp1        
#undef	offset 	    
#undef	a0     	    
#undef	a1     	    
#undef	adr0   	    
#undef	adr1   	    
#undef	brt0  		
#undef	brt1  		
#undef	brt2  		
#undef	brt3  		
#undef	lo0    	    	
#undef	lo1    	    	
#undef	cst_q  	    	
#undef	c_tmp	    	
#undef	a_t0  
#undef	a_t1
#undef	b0     	    	
#undef	b1     	    	
#undef	bp0    	    	
#undef	bp1    	    	
#undef	z0     	    	
#undef	z1     	    	
#undef	i_cnt  	    	
#undef	o_cnt

#undef	Xlo  	    
#undef	Xhi  	    
#undef	Ylo  	    
#undef	Yhi  	    
#undef	Zlo  	    
#undef	Zhi  	    

;epilog
    pop     r31
    pop     r30
    pop     r29
    pop     r28
    pop     r27
    pop     r26
    pop     r25
    pop     r24
    pop     r23
    pop     r22
    pop     r21
    pop     r20
    pop     r19
    pop     r18
    pop     r17
    pop     r16
    pop     r15
    pop     r14
    pop     r13
    pop     r12
    pop     r11
    pop     r10
    pop     r9
    pop     r8
    pop     r7
    pop     r6
    pop     r5
    pop     r4
    pop     r3
    pop     r2
    pop     r1
    pop     r0
    
    ret

    END
    

