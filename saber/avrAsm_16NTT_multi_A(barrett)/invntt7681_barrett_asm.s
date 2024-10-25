
/*!
 * @details     void q2_invntt_mont_asm (int16_t a[256]);
 * @brief       inverse NTT 7681 
 * @param[in]   signed 16-bit x 256 polynomial data : a
 * @param[out]  return invNTT(a)
 * @note        ..
 * @pre         ..
 * @author Youngbeom Kim
 */
            EXTERN  q2_streamlined_CT_asm_inv_table
            EXTERN  q2_last_twist_table
            NAME    q2_barrett_invntt_asm   ;   module name
            PUBLIC  q2_barrett_invntt_asm   ;   make the main label vissible outside this module
            RSEG    CODE                    ;   place program in 'CODE' segment

mc_muls16x16_32    MACRO    m_a0, m_a1, m_b0, m_b1, m_c0, m_c1, m_c2, m_c3, m_zero

	muls    m_a1,   m_b1        ;   (signed)ah * (signed)
	movw	m_c2,   r0
	mul     m_a0,   m_b0        ;   al * bl
	movw	m_c0,   r0
	mulsu	m_a1,   m_b0        ;   (signed)ah * bl
	sbc     m_c3,   m_zero
	add     m_c1,   r0
	adc     m_c2,   r1
	adc     m_c3,   m_zero
	mulsu	m_b1,   m_a0        ;   (signed)bh * al
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
    
mc_load_zetas   MACRO

    lpm     z0,      Z+
    lpm     z1,      Z+
    lpm     bp0,     Z+
    lpm     bp1,     Z+
    
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
    
mc_ct_brt_lo    MACRO    m_a0, m_a1, m_b0, m_b1, m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero
    
    movw    m_s1lo, m_a0
    add     m_a0,   m_b0
    adc     m_a1,   m_b1
    sub     m_s1lo, m_b0
    sbc     m_s1hi, m_b1
    movw    m_b0,   m_s1lo
    
    // a[j] = q2_barrett_reduce((a[j] + t));
    movw    m_s1lo, m_a0
    ldi     m_s2lo, 0x21
    ldi     m_s2hi, 0x22    ; m_s2 = v = 8737 = 0x2221
    mc_muls16x16_32 m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero
    ldi     m_s1lo, 2
    add     m_mont3, m_s1lo
    mov     m_mont2, m_mont3
    lsl     m_mont3
    sbc     m_mont3, m_mont3
    asr     m_mont2
    asr     m_mont2

    movw    m_s1lo, m_mont2
    ldi     m_s2lo, 0x01
    ldi     m_s2hi, 0x1E    ; m_s2 = q = 7681 = 0x1E01
    mc_muls16x16_32 m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero
    sub     m_a0,   m_mont0
    sbc     m_a1,   m_mont1
    st      X,      m_a1
    st      -X,     m_a0
                
    // a[j + len] = q2_barrett_reduce(a[j] - t);
    movw    m_s1lo, m_b0
    ldi     m_s2lo, 0x21
    ldi     m_s2hi, 0x22    ; m_s2 = v = 8737 = 0x2221
    mc_muls16x16_32 m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero   
    ldi     m_s1lo, 2
    add     m_mont3, m_s1lo
    mov     m_mont2, m_mont3
    lsl     m_mont3
    sbc     m_mont3, m_mont3
    asr     m_mont2
    asr     m_mont2

    movw    m_s1lo, m_mont2
    ldi     m_s2lo, 0x01
    ldi     m_s2hi, 0x1E    ; m_s2 = q = 7681 = 0x1E01
    mc_muls16x16_32 m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero
    sub     m_b0,   m_mont0
    sbc     m_b1,   m_mont1

    st      Y+,     m_b0
    st      Y+,     m_b1
    adiw    r26,    2
   
    ENDM    
                      
q2_barrett_invntt_asm:
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
    movw    adr0,   r16
    
    ;   Z = zetas || zetashi table
    ldi     r30,    low (q2_streamlined_CT_asm_inv_table)
    ldi     r31,    high(q2_streamlined_CT_asm_inv_table)
    
layer7: 
    ldi     o_cnt,  128      ;    
    movw    Xlo,    adr0
    movw    Ylo,    adr0
    ldi     c_tmp,  2        ;
    mov     offset, c_tmp
    add     Ylo,    c_tmp
    adc     Yhi,    zero
        
    
out7:            
    mc_ct_hi    a0, a1, b0, b1       
    mc_ct_lo    a0, a1, b0, b1, brt0, brt1    
    
    //end
    
    add     Xlo,    offset
    adc     Xhi,    zero
    add     Ylo,    offset
    adc     Yhi,    zero
    dec     o_cnt
    brne    out7


layer6: // Layer 6    
    ldi     o_cnt,  64       ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  4        ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero           
        
in6:
    mc_ct_hi        a0, a1, b0, b1       
    mc_ct_lo        a0, a1, b0, b1, brt0, brt1 
    
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    sbiw    Zlo,    4
    
    add     Xlo,    offset
    adc     Xhi,    zero
    add     Ylo,    offset
    adc     Yhi,    zero
    
    dec     o_cnt
    breq    endin6          ;
    rjmp    in6
    
endin6:    
    
layer5:
    ldi     o_cnt,  32       ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  8        ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero   
    
    adiw    Zlo,    4        ;
             
in5:
    mc_ct_hi        a0, a1, b0, b1       
    mc_ct_brt_lo    a0, a1, b0, b1, cst_q, c_tmp, bp0, bp1, brt0, brt1, brt2, brt3, zero
    
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero            
    mc_ct_brt_lo    a0, a1, lo0, lo1, b0, b1, bp0, bp1, brt0, brt1, brt2, brt3, zero
        
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_brt_lo    a0, a1, lo0, lo1, b0, b1, bp0, bp1, brt0, brt1, brt2, brt3, zero
        
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_brt_lo    a0, a1, lo0, lo1, b0, b1, bp0, bp1, brt0, brt1, brt2, brt3, zero
    sbiw    Zlo,    12    
    
    add     Xlo,    offset
    adc     Xhi,    zero
    add     Ylo,    offset
    adc     Yhi,    zero
    
    dec     o_cnt
    breq    endin5          ;
    rjmp    in5             ;
    
endin5:

layer4:
    ldi     o_cnt,  16       ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  16        ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero   
    
    adiw    Zlo,    12        ;
             
in4:
    mc_ct_hi        a0, a1, b0, b1       
    mc_ct_lo        a0, a1, b0, b1, brt0, brt1 
    
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
        
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
        
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
        
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
        
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    sbiw    Zlo,    28    
    
    add     Xlo,    offset
    adc     Xhi,    zero
    add     Ylo,    offset
    adc     Yhi,    zero
    
    dec     o_cnt
    breq    endin4          ;
    rjmp    in4             ;
    
endin4:   

layer3:
    ldi     o_cnt,  8        ;
    ldi     i_cnt,  15       ;
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  32        ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero   
    
    adiw    Zlo,    28        ;
             
out3:
    mc_ct_hi        a0, a1, b0, b1       
    mc_ct_lo        a0, a1, b0, b1, brt0, brt1 

in3:
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    dec     i_cnt
    breq    endin3          ;
    rjmp    in3             ;
    
endin3:    
    ldi     i_cnt,  15      ;
    sbiw    Zlo,    60      ;     
    add     Xlo,    offset
    adc     Xhi,    zero
    add     Ylo,    offset
    adc     Yhi,    zero
    
    dec     o_cnt
    breq    endout3          ;
    rjmp    out3             ;
    
endout3:


layer2:
    ldi     o_cnt,  4        ;
    ldi     i_cnt,  31       ;
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  64        ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero   
    
    adiw    Zlo,    60        ;
             
out2:
    mc_ct_hi        a0, a1, b0, b1       
    mc_ct_lo        a0, a1, b0, b1, brt0, brt1 

in2:
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero        
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    dec     i_cnt
    breq    endin2          ;
    rjmp    in2             ;
    
endin2:    
    ldi     i_cnt,  31      ;
    ldi     c_tmp,  124     ;
    sub     Zlo,    c_tmp
    sbc     Zhi,    zero
    add     Xlo,    offset
    adc     Xhi,    zero
    add     Ylo,    offset
    adc     Yhi,    zero
    
    dec     o_cnt
    breq    endout2          ;
    rjmp    out2             ;
    
endout2:

layer1:
    ldi     o_cnt,  2        ;
    ldi     i_cnt,  63       ;
    movw    Xlo,    adr0
    movw    Ylo,    adr0
    ldi     c_tmp,  128      ;
    mov     offset, c_tmp
    add     Ylo,    c_tmp
    adc     Yhi,    zero   
    
    ldi     c_tmp,  124      ;
    add     Zlo,    c_tmp 
    adc     Zhi,    zero  
             
out1:
    mc_ct_hi        a0, a1, b0, b1       
    mc_ct_brt_lo    a0, a1, b0, b1, cst_q, c_tmp, bp0, bp1, brt0, brt1, brt2, brt3, zero

in1:
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero
    mc_ct_brt_lo    a0, a1, lo0, lo1, b0, b1, bp0, bp1, brt0, brt1, brt2, brt3, zero    
    dec     i_cnt
    breq    endin1          ;
    rjmp    in1             ;
    
endin1:    
    ldi     i_cnt,  63      ;
    ldi     cst_q,  252     ;
    sub     Zlo,    cst_q    
    sbc     Zhi,    zero
    add     Xlo,    offset  
    adc     Xhi,    zero    
    add     Ylo,    offset  
    adc     Yhi,    zero
    
    dec     o_cnt
    breq    endout1          ;
    rjmp    out1             ;
    
endout1:


layer0:    
    ldi     i_cnt,  127       ;
    movw    Xlo,    adr0
    movw    Ylo,    adr0
    inc     Yhi
              
    ldi     cst_q,  252     ;
    add     Zlo,    cst_q    
    adc     Zhi,    zero  
             
    mc_ct_hi        a0, a1, b0, b1       
    mc_ct_lo        a0, a1, b0, b1, brt0, brt1 

in0:
    mc_load_zetas
    mc_ct_hi        a0, a1, b0, b1   
    mc_barrett_mul  b0, b1, z0, z1, bp0, bp1, lo0, lo1, brt0, brt1, brt2, brt3, c_tmp, cst_q, zero
    mc_ct_lo        a0, a1, lo0, lo1, brt0, brt1
    dec     i_cnt
    breq    endin0          ;
    rjmp    in0             ;
    
endin0:       
    

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
    