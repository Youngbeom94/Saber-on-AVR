
/*!
 * @details     void q1_invntt_mont_asm (int16_t a[256]);
 * @brief       inverse NTT 7681 
 * @param[in]   signed 16-bit x 256 polynomial data : a
 * @param[out]  return invNTT(a)
 * @note        ..
 * @pre         ..
 * @author Youngbeom Kim
 */
            EXTERN  q1_streamlined_CT_inv_table
            EXTERN  q1_last_twist_mont_table
            NAME    q1_invntt_mont_asm      ;   module name
            PUBLIC  q1_invntt_mont_asm      ;   make the main label vissible outside this module
            RSEG    CODE                 ;   place program in 'CODE' segment

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
    
mc_ct_hi        MACRO    m_a0, m_a1, m_b0, m_b1

	ld      m_a0,   X+
    ld      m_a1,   X
    ldd     m_b0,   Y+0
    ldd     m_b1,   Y+1

    ENDM
    
    
mc_montmul      MACRO    m_b0, m_b1, m_z0, m_z1, m_lo0, m_lo1, m_hi0, m_hi1, m_mont0, m_mont1, m_mont2, m_mont3, m_cst_tmp, m_cst_q, m_zero

    ldi     m_cst_q,    0xE2 ; QinvR
    mc_muls16x16_32 m_b0, m_b1, m_z0, m_z1, m_lo0, m_lo1, m_hi0, m_hi1, m_zero
    mulsu   m_cst_q,    m_lo0
    add     m_lo1,      r0
    
    ldi     m_cst_q,    0x1E ; Qhi
    ldi     m_cst_tmp,  0x01 ; Qlo
    mc_muls16x16_32 m_lo0, m_lo1, m_cst_tmp, m_cst_q, m_mont0, m_mont1, m_mont2, m_mont3, m_zero
    sub     m_hi0,      m_mont2
    sbc     m_hi1,      m_mont3
    
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
    
mc_ct_barrett_lo        MACRO    m_a0, m_a1, m_b0, m_b1, m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero
    
    movw    m_s1lo, m_a0
    add     m_a0,   m_b0
    adc     m_a1,   m_b1
    sub     m_s1lo, m_b0
    sbc     m_s1hi, m_b1
    movw    m_b0,   m_s1lo
    
    // a[j] = q1_barrett_reduce((a[j] + t));
    movw    m_s1lo, m_a0
    ldi     m_s2lo, 0x21
    ldi     m_s2hi, 0x22    ; m_s2 = v = 8737 = 0x2221
    mc_muls16x16_32 m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero   
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
                
    // a[j + len] = q1_barrett_reduce(a[j] - t);
    movw    m_s1lo, m_b0
    ldi     m_s2lo, 0x21
    ldi     m_s2hi, 0x22    ; m_s2 = v = 8737 = 0x2221
    mc_muls16x16_32 m_s1lo, m_s1hi, m_s2lo, m_s2hi, m_mont0, m_mont1, m_mont2, m_mont3, m_zero   
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
                      
q1_invntt_mont_asm:
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

#define	zero    r2
#define	offset 	r5	
#define	a0     	r6	
#define	a1     	r7	
#define	adr0   	r8	
#define	adr1   	r9	
#define	mont0  	r10	
#define	mont1  	r11	
#define	mont2  	r12	
#define	mont3  	r13	
#define	hi0    	r14	
#define	hi1    	r15	
#define	cst_q  	r16	
#define	c_tmp	r17	
#define	b0     	r18	
#define	b1     	r19	
#define	lo0    	r20	
#define	lo1    	r21	
#define	z0     	r22	
#define	z1     	r23	
#define	i_cnt  	r24	
#define	o_cnt  	r25	

initialize:
    clr     zero
    ;   X = (int16_t*)&a, Y= (int16_t*)&a + 128
    movw    adr0,   r16
    
    ;   Z = zetas || zetashi table
    ldi     r30,    low (q1_streamlined_CT_inv_table)
    ldi     r31,    high(q1_streamlined_CT_inv_table)
    
layer7: 
    ldi     o_cnt,  128      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  2        ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
    
    lpm     z0,     Z+
    lpm     z1,     Z+
    
out7:        
    
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_lo    a0, a1, hi0, hi1, mont0, mont1    
    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    dec     o_cnt
    brne    out7


layer6: // Layer 6
    ldi     i_cnt,  2       ;
    ldi     o_cnt,  64      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  4       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
             
in6:
    lpm     z0,     Z+
    lpm     z1,     Z+
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_lo    a0, a1, hi0, hi1, mont0, mont1           
    
    dec     i_cnt
    breq    endin6          ;
    rjmp    in6
    
endin6:    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    
    ldi     c_tmp,  4       ;
    sub     r30,    c_tmp
    sbc     r31,    zero
    
    ldi     i_cnt,  2      ;
    
    dec     o_cnt
    breq    layer5  ;
    rjmp    in6
    
layer5:
    ldi     c_tmp,  4       ;
    add     r30,    c_tmp
    adc     r31,    zero
    
    ldi     i_cnt,  4       ;
    ldi     o_cnt,  32      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  8       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
             
in5:
    lpm     z0,     Z+
    lpm     z1,     Z+
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_lo    a0, a1, hi0, hi1, mont0, mont1           
    
    dec     i_cnt
    breq    endin5          ;
    rjmp    in5             ;
    
endin5:    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    
    ldi     c_tmp,  8       ;
    sub     r30,    c_tmp
    sbc     r31,    zero
    
    ldi     i_cnt,  4       ;
    
    dec     o_cnt
    breq    layer4          ;
    rjmp    in5             ;

layer4:
    ldi     c_tmp,  8       ;
    add     r30,    c_tmp
    adc     r31,    zero
    
    ldi     i_cnt,  8       ;
    ldi     o_cnt,  16      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  16       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
             
in4:
    lpm     z0,     Z+
    lpm     z1,     Z+
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_lo    a0, a1, hi0, hi1, mont0, mont1           
    
    dec     i_cnt
    breq    endin4          ;
    rjmp    in4             ;
    
endin4:    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    
    ldi     c_tmp,  16       ;
    sub     r30,    c_tmp
    sbc     r31,    zero
    
    ldi     i_cnt,  8       ;
    
    dec     o_cnt
    breq    layer3          ;
    rjmp    in4             ;

layer3:
    ldi     c_tmp,  16       ;
    add     r30,    c_tmp
    adc     r31,    zero
    
    ldi     i_cnt,  16       ;
    ldi     o_cnt,  8      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  32       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
             
in3:
    lpm     z0,     Z+
    lpm     z1,     Z+
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_lo    a0, a1, hi0, hi1, mont0, mont1           
    
    dec     i_cnt
    breq    endin3          ;
    rjmp    in3             ;
    
endin3:    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    
    ldi     c_tmp,  32       ;
    sub     r30,    c_tmp
    sbc     r31,    zero
    
    ldi     i_cnt,  16       ;
    
    dec     o_cnt
    breq    layer2          ;
    rjmp    in3             ;

layer2:
    ldi     c_tmp,  32       ;
    add     r30,    c_tmp
    adc     r31,    zero
    
    ldi     i_cnt,  32       ;
    ldi     o_cnt,  4      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  64       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
             
in2:
    lpm     z0,     Z+
    lpm     z1,     Z+
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_barrett_lo    a0, a1, hi0, hi1, lo0, lo1, cst_q, c_tmp, mont0, mont1, mont2, mont3, zero       
    
    dec     i_cnt
    breq    endin2          ;
    rjmp    in2             ;
    
endin2:    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    
    ldi     c_tmp,  64       ;
    sub     r30,    c_tmp
    sbc     r31,    zero
    
    ldi     i_cnt,  32       ;
    
    dec     o_cnt
    breq    layer1          ;
    rjmp    in2             ;

layer1:
    ldi     c_tmp,  64       ;
    add     r30,    c_tmp
    adc     r31,    zero
    
    ldi     i_cnt,  64       ;
    ldi     o_cnt,  2      ;    
    movw    r26,    adr0
    movw    r28,    adr0
    ldi     c_tmp,  128       ;
    mov     offset, c_tmp
    add     r28,    c_tmp
    adc     r29,    zero
             
in1:
    lpm     z0,     Z+
    lpm     z1,     Z+
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_lo    a0, a1, hi0, hi1, mont0, mont1           
    
    dec     i_cnt
    breq    endin1          ;
    rjmp    in1             ;
    
endin1:    
    add     r26,    offset
    adc     r27,    zero
    add     r28,    offset
    adc     r29,    zero
    
    ldi     c_tmp,  128      ;
    sub     r30,    c_tmp
    sbc     r31,    zero
    
    ldi     i_cnt,  64       ;
    
    dec     o_cnt
    breq    layer0          ;
    rjmp    in1             ;

layer0:
    ldi     c_tmp,  128       ;
    add     r30,    c_tmp
    adc     r31,    zero
    
    ldi     i_cnt,  128       ;    
    movw    r26,    adr0
    movw    r28,    adr0    
    inc     r29
             
in0:
    lpm     z0,     Z+
    lpm     z1,     Z+
    mc_ct_hi    a0, a1, b0, b1   
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero
    mc_ct_lo    a0, a1, hi0, hi1, mont0, mont1           
    
    dec     i_cnt
    breq    twist            ;
    rjmp    in0              ;
    
twist:

    ;   Z = zetas || zetashi table
    ldi     r30,    low (q1_last_twist_mont_table)
    ldi     r31,    high(q1_last_twist_mont_table)
    ldi     i_cnt,  0xFF
    
    movw    r28,    adr0
    
twist_in:    
    ldd     b0,     Y+0
    ldd     b1,     Y+1
    lpm     z0,     Z+
    lpm     z1,     Z+
    
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero    
        
    st      Y+,     hi0   
    st      Y+,     hi1
    dec     i_cnt
    brne    twist_in
    
    ldd     b0,     Y+0
    ldd     b1,     Y+1
    lpm     z0,     Z+
    lpm     z1,     Z+
    
    mc_montmul  b0, b1, z0, z1, lo0, lo1, hi0, hi1, mont0, mont1, mont2, mont3, c_tmp, cst_q, zero    
        
    st      Y+,     hi0   
    st      Y+,     hi1
    
end_ntt:

#undef	zero  
#undef	offset
#undef	a0    
#undef	a1    
#undef	adr0  
#undef	adr1  
#undef	mont0 
#undef	mont1 
#undef	mont2 
#undef	mont3 
#undef	hi0   
#undef	hi1   
#undef	cst_q 
#undef	c_tmp	
#undef	b0    
#undef	b1    
#undef	lo0   
#undef	lo1   
#undef	z0    
#undef	z1    
#undef	i_cnt 
#undef	o_cnt 

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
    

