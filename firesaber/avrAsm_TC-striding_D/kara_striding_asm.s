
            EXTERN  res_p00
            EXTERN  res_p01
            EXTERN  res_p02
            EXTERN  res_p10
            EXTERN  res_p11
            EXTERN  res_p12
            EXTERN  res_p20
            EXTERN  res_p21
            EXTERN  res_p22
            
            NAME    karatsuba_striding_asm  ;   module name
            PUBLIC  karatsuba_striding_asm  ;   make the main label vissible outside this module
            RSEG    CODE                    ;   place program in 'CODE' segment
            

mc_add16_res    MACRO    m_rlo, m_rhi, m_s1lo, m_s1hi, m_s2lo, m_s2hi

	movw    m_rlo,  m_s1lo
    add     m_rlo,  m_s2lo
    adc     m_rhi,  m_s2hi

    ENDM
    
mc_mul16x16_lo_acc  MACRO    m_rlo, m_rhi, m_s1lo, m_s1hi, m_s2lo, m_s2hi

	mul     m_s1lo, m_s2lo
    add     m_rlo,  r0
    adc     m_rhi,  r1
    mul     m_s1lo, m_s2hi
    add     m_rhi,  r0
    mul     m_s1hi, m_s2lo
    add     m_rhi,  r0
    
    ENDM
    
mc_mul16x16_lo_scc  MACRO    m_rlo, m_rhi, m_s1lo, m_s1hi, m_s2lo, m_s2hi

	mul     m_s1lo, m_s2lo
    sub     m_rlo,  r0
    sbc     m_rhi,  r1
    mul     m_s1lo, m_s2hi
    sub     m_rhi,  r0
    mul     m_s1hi, m_s2lo
    sub     m_rhi,  r0
    
    ENDM    
    
mc_array_acc  MACRO    m_rlo, m_rhi, m_ofs, m_zero, m_plo, m_phi, m_slo, m_shi, m_solo, m_sohi m_array
    
    ldi     Zlo,    low (m_array)
    ldi     Zhi,    high(m_array)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero
    
    ldd     m_rlo,  Z+0
    ldd     m_rhi,  Z+1

    ldd     m_slo,  Y+m_solo 
    ldd     m_shi,  Y+m_sohi
    
    mc_mul16x16_lo_acc  m_rlo, m_rhi, m_plo, m_phi, m_slo, m_shi 
    
    std     Z+0,    m_rlo
    std     Z+1,    m_rhi
    
    ENDM
    
mc_light_array_acc  MACRO    m_rlo, m_rhi, m_ofs, m_zero, m_plo, m_phi, m_slo, m_shi, m_array
    
    ldi     Zlo,    low (m_array)
    ldi     Zhi,    high(m_array)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero
    
    ldd     m_rlo,  Z+0
    ldd     m_rhi,  Z+1
    
    mc_mul16x16_lo_acc  m_rlo, m_rhi, m_plo, m_phi, m_slo, m_shi 
    
    std     Z+0,    m_rlo
    std     Z+1,    m_rhi
        
    ENDM
    
mc_array_scc  MACRO    m_rlo, m_rhi, m_ofs, m_zero, m_plo, m_phi, m_slo, m_shi, m_solo, m_sohi m_array
    
    ldi     Zlo,    low (m_array)
    ldi     Zhi,    high(m_array)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero
    
    ldd     m_rlo,  Z+0
    ldd     m_rhi,  Z+1

    ldd     m_slo,  Y+m_solo 
    ldd     m_shi,  Y+m_sohi
    
    mc_mul16x16_lo_scc  m_rlo, m_rhi, m_plo, m_phi, m_slo, m_shi 
    
    std     Z+0,    m_rlo
    std     Z+1,    m_rhi
    
    ENDM

mc_light_array_scc  MACRO    m_rlo, m_rhi, m_ofs, m_zero, m_plo, m_phi, m_slo, m_shi, m_array
    
    ldi     Zlo,    low (m_array)
    ldi     Zhi,    high(m_array)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero
    
    ldd     m_rlo,  Z+0
    ldd     m_rhi,  Z+1
    
    mc_mul16x16_lo_scc  m_rlo, m_rhi, m_plo, m_phi, m_slo, m_shi 
    
    std     Z+0,    m_rlo
    std     Z+1,    m_rhi
        
    ENDM    

mc_3ld_scc  MACRO   m_reslo, m_reshi, m_resary, m_s1lo, m_s1hi, m_s1ary, m_s2lo, m_s2hi, m_s2ary, m_ofs, m_zero

    ldi     Zlo,    low (m_s2ary)
    ldi     Zhi,    high(m_s2ary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero    
    ldd     m_s2lo, Z+0
    ldd     m_s2hi, Z+1    
    
    ldi     Zlo,    low (m_s1ary)
    ldi     Zhi,    high(m_s1ary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero    
    ldd     m_s1lo, Z+0
    ldd     m_s1hi, Z+1
    
    ldi     Zlo,    low (m_resary)
    ldi     Zhi,    high(m_resary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero    
    ldd     m_reslo,Z+0
    ldd     m_reshi,Z+1
        
    sub     m_reslo,m_s2lo
    sbc     m_reshi,m_s2hi
    sub     m_reslo,m_s1lo
    sbc     m_reshi,m_s1hi
    
    std     Z+0,    m_reslo
    std     Z+1,    m_reshi
    
    ENDM
    
mc_1ld_acc   MACRO   m_reslo, m_reshi, m_resary, m_s1lo, m_s1hi, m_ofs, m_zero
   
    ldi     Zlo,    low (m_resary)
    ldi     Zhi,    high(m_resary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero    
    ldd     m_reslo,Z+0
    ldd     m_reshi,Z+1
        
    add     m_reslo,m_s1lo
    adc     m_reshi,m_s1hi
    
    std     Z+0,    m_reslo
    std     Z+1,    m_reshi
    
    ENDM

mc_2ld_sub_Xst  MACRO   m_s1lo, m_s1hi, m_s1ary, m_s2lo, m_s2hi, m_s2ary, m_ofs, m_zero
   
    ldi     Zlo,    low (m_s1ary)
    ldi     Zhi,    high(m_s1ary)
    ldd     m_s1lo, Z+0
    ldd     m_s1hi, Z+1
    
    ldi     Zlo,    low (m_s2ary)
    ldi     Zhi,    high(m_s2ary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero    
    ldd     m_s2lo, Z+0
    ldd     m_s2hi, Z+1
                
    sub     m_s1lo, m_s2lo
    sbc     m_s1hi, m_s2hi
    
    st      X+,     m_s1lo
    st      X+,     m_s1hi
    
    ENDM

mc_1ld_Xst  MACRO   m_s1lo, m_s1hi, m_s1ary, m_ofs, m_zero
           
    ldi     Zlo,    low (m_s1ary)
    ldi     Zhi,    high(m_s1ary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero    
    ldd     m_s1lo, Z+0
    ldd     m_s1hi, Z+1          
    
    st      X+,     m_s1lo
    st      X+,     m_s1hi
    
    ENDM
    
mc_1ld_bfst MACRO   m_bflo, m_bfhi, m_s1lo, m_s1hi, m_s1ary, m_ofs, m_zero
           
    ldi     Zlo,    low (m_s1ary)
    ldi     Zhi,    high(m_s1ary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero    
    ldd     m_s1lo, Z+0
    ldd     m_s1hi, Z+1          
    
    movw    m_bflo, m_s1lo
    
    ENDM
    
    
mc_1ld1bf_add_Xst   MACRO   m_s1lo, m_s1hi, m_s1ary, m_bflo, m_bfhi, m_ofs, m_zero
   
    ldi     Zlo,    low (m_s1ary)
    ldi     Zhi,    high(m_s1ary)
    add     Zlo,    m_ofs
    adc     Zhi,    m_zero 
    ldd     m_s1lo, Z+0
    ldd     m_s1hi, Z+1
                   
    add     m_bflo, m_s1lo
    adc     m_bfhi, m_s1hi
    
    st      X+,     m_bflo
    st      X+,     m_bfhi
    
    ENDM
    
karatsuba_striding_asm:

//  prolog
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

#define	p00_lo  r2  
#define	p00_hi 	r3  
#define	p01_lo 	r4  
#define	p01_hi 	r5  
#define	p02_lo 	r6  
#define	p02_hi 	r7  
#define	p10_lo 	r8	
#define	p10_hi 	r9	
#define	p11_lo 	r10	
#define	p11_hi 	r11	
#define	p12_lo 	r12	
#define	p12_hi 	r13	
#define	p20_lo 	r14	
#define	p20_hi	r15	
#define	p21_lo 	r16	
#define	p21_hi 	r17	
#define	p22_lo 	r18	
#define	p22_hi 	r19
#define	tp0_lo 	r20
#define	tp0_hi 	r21
#define	tp3_lo 	r22
#define	tp3_hi 	r23
#define	tp2_lo  r24
#define	tp2_hi  r25
#define	tp1_lo 	r26
#define	tp1_hi 	r27

#define	loop_i  r24
#define	loop_o  r25

#define	Xlo 	r26
#define	Xhi 	r27
#define	Ylo 	r28
#define	Yhi 	r29
#define	Zlo 	r30
#define	Zhi 	r31

#define bf0_lo  r0
#define bf0_hi  r1
#define bf1_lo  r20
#define bf1_hi  r21
#define bf2_lo  r28
#define bf2_hi  r29

initialize:
    
    movw    Xlo,    r16
    movw    Ylo,    r20
    push    Xlo                                                     ; src1_deslo -> stack 1
    push    Xhi                                                     ; src1_deshi -> stack 2
    push    Xlo                                                     ; src1_deslo -> stack 3
    push    Xhi                                                     ; src1_deshi -> stack 4

    clr     loop_o
    
    //  Load 4 and Cache 9
LOOP_O1:
    pop     Xhi                                                     ; stack -> src1_deshi 3
    pop     Xlo                                                     ; stack -> src1_deslo 2
    
    ld      p00_lo, X+      ; p00 = src1_des[4 * i + 0];
    ld      p00_hi, X+
    
    ld      p02_lo, X+      ; p02 = src1_des[4 * i + 1];
    ld      p02_hi, X+
    
    ld      p20_lo, X+      ; p20 = src1_des[4 * i + 2];
    ld      p20_hi, X+
    
    ld      p22_lo, X+      ; p22 = src1_des[4 * i + 3];
    ld      p22_hi, X+
    
    push    Xlo                                                     ; src1_deslo -> stack 3
    push    Xhi                                                     ; src1_deshi -> stack 4

    mc_add16_res    p01_lo, p01_hi, p00_lo, p00_hi, p02_lo, p02_hi  ; p01 = p00 + p02;
    mc_add16_res    p21_lo, p21_hi, p20_lo, p20_hi, p22_lo, p22_hi  ; p21 = p20 + p22;
    
    mc_add16_res    p10_lo, p10_hi, p00_lo, p00_hi, p20_lo, p20_hi  ; p10 = p00 + p20;
    mc_add16_res    p12_lo, p12_hi, p02_lo, p02_hi, p22_lo, p22_hi  ; p12 = p02 + p22;
    
    mc_add16_res    p11_lo, p11_hi, p10_lo, p10_hi, p12_lo, p12_hi  ; p11 = p10 + p12;
    
    
    push    loop_o                                                  ; i -> stack 5
    
    mov     tp1_lo, loop_o  ; tp1_lo = i
    
    ldi     loop_i, 16
    sub     loop_i, tp1_lo  ; loop_i = 16 - i
    clr     loop_o          ; loop_o = j    

    push    Ylo                                                     ; src2lo -> stack 6
    push    Yhi                                                     ; src2hi -> stack 7
    
    add     tp1_lo,  loop_o ; tp1_lo = i + j
    lsl     tp1_lo
    
LOOP_I1:
    clr     tp1_hi          ; tp1_hi is zero    
    push    loop_i                                                  ; 16 - i -> stack,      tp2_lo On 8
    push    loop_o                                                  ; j -> stack,           tp2_hi On 9
        
    ; res_p00[i + j] += p00 * src2[4 * j + 0], tp2 = src2[4 * j + 0]
    mc_array_acc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p00_lo, p00_hi, tp2_lo, tp2_hi, 0, 1, res_p00
    
    ; res_p02[i + j] += p02 * src2[4 * j + 1], tp3 = src2[4 * j + 1]
    mc_array_acc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p02_lo, p02_hi, tp3_lo, tp3_hi, 2, 3, res_p02
 
    ; _p01  = src2[4 * j + 0] + src2[4 * j + 1], tp3 = _p01
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi
    push    tp3_lo                                                  ; _p01lo -> stack 10
    push    tp3_hi                                                  ; _p01hi -> stack 11
    
    ; res_p01[i + j] += p01 * _p01;,
    mc_light_array_acc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p01_lo, p01_hi, tp3_lo, tp3_hi, res_p01
    
    
    ; res_p20[i + j] += p20 * src2[4 * j + 2], tp2 = src2[4 * j + 2]
    mc_array_acc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p20_lo, p20_hi, tp2_lo, tp2_hi, 4, 5, res_p20
    
    ; res_p22[i + j] += p22 * src2[4 * j + 3], tp3 = src2[4 * j + 3]
    mc_array_acc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p22_lo, p22_hi, tp3_lo, tp3_hi, 6, 7, res_p22
 
    ; _p21  = src2[4 * j + 2] + src2[4 * j + 3], tp3 = _p01
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi
    push    tp3_lo                                                  ; _p21lo -> stack 12
    push    tp3_hi                                                  ; _p21hi -> stack 13
    
    ; res_p01[i + j] += p21 * _p21;,
    mc_light_array_acc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p21_lo, p21_hi, tp3_lo, tp3_hi, res_p21
    
    
    ; res_p10[i + j] += p10 * (src2[4 * j + 0] + src2[4 * j + 2]);
    ldd     tp2_lo, Y+0     ; tp2 =  src2[4 * j + 0]
    ldd     tp2_hi, Y+1
    ldd     tp3_lo, Y+4     ; tp3 =  src2[4 * j + 2]
    ldd     tp3_hi, Y+5
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi  ; tp3 = src2[4 * j + 0] + src2[4 * j + 2]       
    mc_light_array_acc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p10_lo, p10_hi, tp3_lo, tp3_hi, res_p10
    
    ; res_p12[i + j] += p12 * (src2[4 * j + 1] + src2[4 * j + 3]);
    ldd     tp2_lo, Y+2     ; tp2 =  src2[4 * j + 1]
    ldd     tp2_hi, Y+3
    ldd     tp3_lo, Y+6     ; tp3 =  src2[4 * j + 3]
    ldd     tp3_hi, Y+7
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi  ; tp3 = src2[4 * j + 1] + src2[4 * j + 3]       
    mc_light_array_acc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p12_lo, p12_hi, tp3_lo, tp3_hi, res_p12
    
    ; _p11 = _p01 + _p21
    pop     tp3_hi                                                  ; stack -> _p21hi 12
    pop     tp3_lo                                                  ; stack -> _p21lo 11
    pop     tp2_hi                                                  ; stack -> _p01hi 10
    pop     tp2_lo                                                  ; stack -> _p01lo  9
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi
    
    ; res_p11[i + j] += p11 * _p11;
    mc_light_array_acc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p11_lo, p11_hi, tp3_lo, tp3_hi, res_p11
    
    adiw    Ylo,    8       ; Y + 4 * 2 // control offset 
    
    pop     loop_o                                                  ; stack -> j,           tp2_lo Off 8
    pop     loop_i                                                  ; stack -> 16 - i,      tp2_lo Off 7
    inc     loop_o
    inc     tp1_lo          ; tp1_lo = i + j
    inc     tp1_lo      
    cp      loop_o, loop_i    
    breq    LOOP_I1END
    rjmp    LOOP_I1
    
LOOP_I1END:

    pop     Yhi                                                     ; stack -> src2_hi 6
    pop     Ylo                                                     ; stack -> src2_lo 5
    pop     loop_o                                                  ; stack -> i 4
    inc     loop_o          ; loop_o = i + 1
    
    clr     loop_i
    inc     loop_i
    cp      loop_o, loop_i  
    brne    LOOP_I2ENABLE
    rjmp    LOOP_O1

LOOP_I2ENABLE:
    dec     loop_o          ; loop_o = i
    
    push    loop_o                                                  ; stack -> i 5
    push    Ylo                                                     ; src2lo -> stack 6
    push    Yhi                                                     ; src2hi -> stack 7 

    ldi     loop_i, 16
    sub     loop_i, loop_o  ; loop_i = j =  16 - i

    add     loop_o, loop_i  ; Loop_o = i + j
    ldi     tp1_lo, 16
    sub     loop_o, tp1_lo  ; Loop_o = i + j - 16
    mov     tp1_lo, loop_o  ; tp1_lo = [i + j - 16]
    lsl     tp1_lo
    
    clr     tp1_hi
    mov     tp0_lo, loop_i
    lsl     tp0_lo
    lsl     tp0_lo
    lsl     tp0_lo
    
    add     Ylo,    tp0_lo
    adc     Yhi,    tp1_hi        
    
LOOP_I2:
    push    loop_i                                                  ; j -> stack,      tp2_lo On 8
    
    ; res_p00[i + j - 16] -= p00 * src2[4 * j + 0], tp2 = src2[4 * j + 0]
    mc_array_scc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p00_lo, p00_hi, tp2_lo, tp2_hi, 0, 1, res_p00    

    ; res_p02[i + j - 16] -= p02 * src2[4 * j + 1], tp3 = src2[4 * j + 1]
    mc_array_scc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p02_lo, p02_hi, tp3_lo, tp3_hi, 2, 3, res_p02
 
    ; _p01  = src2[4 * j + 0] + src2[4 * j + 1], tp3 = _p01
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi
    push    tp3_lo                                                  ; _p01lo -> stack 10
    push    tp3_hi                                                  ; _p01hi -> stack 11
    
    ; res_p01[i + j - 16] -= p01 * _p01;,
    mc_light_array_scc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p01_lo, p01_hi, tp3_lo, tp3_hi, res_p01

    ; res_p20[i + j- 16] -= p20 * src2[4 * j + 2], tp2 = src2[4 * j + 2]
    mc_array_scc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p20_lo, p20_hi, tp2_lo, tp2_hi, 4, 5, res_p20
    
    ; res_p22[i + j- 16] -= p22 * src2[4 * j + 3], tp3 = src2[4 * j + 3]
    mc_array_scc        tp0_lo, tp0_hi, tp1_lo, tp1_hi, p22_lo, p22_hi, tp3_lo, tp3_hi, 6, 7, res_p22
 
    ; _p21  = src2[4 * j + 2] + src2[4 * j + 3], tp3 = _p01
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi
    push    tp3_lo                                                  ; _p21lo -> stack 12
    push    tp3_hi                                                  ; _p21hi -> stack 13
    
    ; res_p01[i + j- 16] -= p21 * _p21;,
    mc_light_array_scc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p21_lo, p21_hi, tp3_lo, tp3_hi, res_p21
    
    ; res_p10[i + j] -= p10 * (src2[4 * j + 0] + src2[4 * j + 2]);
    ldd     tp2_lo, Y+0     ; tp2 =  src2[4 * j + 0]
    ldd     tp2_hi, Y+1
    ldd     tp3_lo, Y+4     ; tp3 =  src2[4 * j + 2]
    ldd     tp3_hi, Y+5
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi  ; tp3 = src2[4 * j + 0] + src2[4 * j + 2]       
    mc_light_array_scc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p10_lo, p10_hi, tp3_lo, tp3_hi, res_p10
    
    ; res_p12[i + j] -= p12 * (src2[4 * j + 1] + src2[4 * j + 3]);
    ldd     tp2_lo, Y+2     ; tp2 =  src2[4 * j + 1]
    ldd     tp2_hi, Y+3
    ldd     tp3_lo, Y+6     ; tp3 =  src2[4 * j + 3]
    ldd     tp3_hi, Y+7
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi  ; tp3 = src2[4 * j + 1] + src2[4 * j + 3]       
    mc_light_array_scc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p12_lo, p12_hi, tp3_lo, tp3_hi, res_p12
    
    ; _p11 = _p01 + _p21
    pop     tp3_hi                                                  ; stack -> _p21hi 12
    pop     tp3_lo                                                  ; stack -> _p21lo 11
    pop     tp2_hi                                                  ; stack -> _p01hi 10
    pop     tp2_lo                                                  ; stack -> _p01lo  9
    add     tp3_lo, tp2_lo
    adc     tp3_hi, tp2_hi
    
    ; res_p11[i + j] -= p11 * _p11;
    mc_light_array_scc  tp0_lo, tp0_hi, tp1_lo, tp1_hi, p11_lo, p11_hi, tp3_lo, tp3_hi, res_p11

    adiw    Ylo,    8       ; Y + 4 * 2 // control offset 

    inc     tp1_lo
    inc     tp1_lo
    ldi     tp0_lo, 16
    pop     loop_i                                                  ; stack -> j,      tp2_lo Off 7
    inc     loop_i
    cp      loop_i, tp0_lo
    breq    LOOP_I2END
    rjmp    LOOP_I2
    
LOOP_I2END:
    
    pop     Yhi                                                     ; stack -> src2_hi 6
    pop     Ylo                                                     ; stack -> src2_lo 5
    pop     loop_o                                                  ; stack -> i 4
    inc     loop_o
    ldi     tp0_lo, 16  
    cp      loop_o, tp0_lo   
    breq    LOOP_O1END
    rjmp    LOOP_O1
      
LOOP_O1END:
    
    pop     Xhi                                                     ; stack -> src1_deslo 3
    pop     Xlo                                                     ; stack -> src1_deshi 2

    clr     loop_o
    ldi     loop_i, 16    
    
LOOP_O2:
    mov     Ylo,    loop_o  ; Ylo is offset
    lsl     Ylo
    clr     Yhi             ; Yhi is zero
    
    mc_3ld_scc  p01_lo, p01_hi, res_p01, p00_lo, p00_hi, res_p00, p02_lo, p02_hi, res_p02, Ylo, Yhi ; res_p01[i] = res_p01[i] - res_p00[i] - res_p02[i];
    mc_3ld_scc  p11_lo, p11_hi, res_p11, p10_lo, p10_hi, res_p10, p12_lo, p12_hi, res_p12, Ylo, Yhi ; res_p11[i] = res_p11[i] - res_p10[i] - res_p12[i];
    mc_3ld_scc  p21_lo, p21_hi, res_p21, p20_lo, p20_hi, res_p20, p22_lo, p22_hi, res_p22, Ylo, Yhi ; res_p21[i] = res_p21[i] - res_p20[i] - res_p22[i];
    
    mc_3ld_scc  p10_lo, p10_hi, res_p10, p00_lo, p00_hi, res_p00, p20_lo, p20_hi, res_p20, Ylo, Yhi ; res_p10[i] = res_p10[i] - res_p00[i] - res_p20[i];
    mc_3ld_scc  p11_lo, p11_hi, res_p11, p01_lo, p01_hi, res_p01, p21_lo, p21_hi, res_p21, Ylo, Yhi ; res_p11[i] = res_p11[i] - res_p01[i] - res_p21[i];
    mc_3ld_scc  p12_lo, p12_hi, res_p12, p02_lo, p02_hi, res_p02, p22_lo, p22_hi, res_p22, Ylo, Yhi ; res_p12[i] = res_p12[i] - res_p02[i] - res_p22[i];
    
    mc_1ld_acc  p10_lo, p10_hi, res_p10, p02_lo, p02_hi, Ylo, Yhi                                   ; res_p10[i] += res_p02[i];
    mc_1ld_acc  p20_lo, p20_hi, res_p20, p12_lo, p12_hi, Ylo, Yhi                                   ; res_p20[i] += res_p12[i];
    
    inc     loop_o
    cp      loop_o, loop_i
    breq    LOOP_O2END
    rjmp    LOOP_O2
    
LOOP_O2END:

    pop     Xhi                                                     ; stack -> src1_deslo 1
    pop     Xlo                                                     ; stack -> src1_deshi 0
    
    ldi     tp3_lo, 15      ; tp3_lo = offset
    lsl     tp3_lo
    clr     tp3_hi          ; tp3_hi = zero   
   
    mc_2ld_sub_Xst  p00_lo, p00_hi, res_p00, p20_lo, p20_hi, res_p20, tp3_lo, tp3_hi    ; src1_des[0]  = res_p00[0] - res_p20[15];
    mc_2ld_sub_Xst  p01_lo, p01_hi, res_p01, p21_lo, p21_hi, res_p21, tp3_lo, tp3_hi    ; src1_des[1]  = res_p01[0] - res_p21[15];
    mc_2ld_sub_Xst  p10_lo, p10_hi, res_p10, p22_lo, p22_hi, res_p22, tp3_lo, tp3_hi    ; src1_des[2]  = res_p10[0] - res_p22[15];
    mc_1ld_Xst      p11_lo, p11_hi, res_p11, tp3_hi, tp3_hi                             ; src1_des[3]  = res_p11[0];
    mc_1ld_bfst     bf0_lo, bf0_hi, p20_lo, p20_hi, res_p20, tp3_hi, tp3_hi             ; buff[0] = tst_p20[0];
    mc_1ld_bfst     bf1_lo, bf1_hi, p21_lo, p21_hi, res_p21, tp3_hi, tp3_hi             ; buff[1] = tst_p21[0];
    mc_1ld_bfst     bf2_lo, bf2_hi, p22_lo, p22_hi, res_p22, tp3_hi, tp3_hi             ; buff[2] = tst_p22[0];

    ldi     loop_o, 1
    ldi     loop_i, 15
    clr     tp3_hi          ; tp3_hi = zero 

LOOP_O3:
    mov     tp3_lo, loop_o  ; tp3_lo = offset
    lsl     tp3_lo      
   
    mc_1ld1bf_add_Xst   p00_lo, p00_hi, res_p00, bf0_lo, bf0_hi, tp3_lo, tp3_hi         ; src1_des[4 * i + 0] = buff[0] + res_p00[i];
    mc_1ld1bf_add_Xst   p01_lo, p01_hi, res_p01, bf1_lo, bf1_hi, tp3_lo, tp3_hi         ; src1_des[4 * i + 1] = buff[1] + res_p01[i];
    mc_1ld1bf_add_Xst   p10_lo, p10_hi, res_p10, bf2_lo, bf2_hi, tp3_lo, tp3_hi         ; src1_des[4 * i + 2] = buff[2] + res_p10[i];
    mc_1ld_Xst          p11_lo, p11_hi, res_p11, tp3_lo, tp3_hi                         ; src1_des[4 * i + 3] = res_p11[i];
    mc_1ld_bfst         bf0_lo, bf0_hi, p20_lo, p20_hi, res_p20, tp3_lo, tp3_hi         ; buff[0] = res_p20[i];
    mc_1ld_bfst         bf1_lo, bf1_hi, p21_lo, p21_hi, res_p21, tp3_lo, tp3_hi         ; buff[1] = res_p21[i];
    mc_1ld_bfst         bf2_lo, bf2_hi, p22_lo, p22_hi, res_p22, tp3_lo, tp3_hi         ; buff[2] = res_p22[i];

    inc     loop_o
    cp      loop_o, loop_i
    breq    LOOP_O3END
    rjmp    LOOP_O3

LOOP_O3END:

    mov     tp3_lo, loop_o  ; tp3_lo = offset
    lsl     tp3_lo
    
    mc_1ld1bf_add_Xst   p00_lo, p00_hi, res_p00, bf0_lo, bf0_hi, tp3_lo, tp3_hi         ; src1_des[60] = buff[0] + res_p00[15];
    mc_1ld1bf_add_Xst   p01_lo, p01_hi, res_p01, bf1_lo, bf1_hi, tp3_lo, tp3_hi         ; src1_des[61] = buff[1] + res_p01[15];
    mc_1ld1bf_add_Xst   p10_lo, p10_hi, res_p10, bf2_lo, bf2_hi, tp3_lo, tp3_hi         ; src1_des[62] = buff[2] + res_p10[15];
    mc_1ld_Xst          p11_lo, p11_hi, res_p11, tp3_lo, tp3_hi                         ; src1_des[63] = res_p11[15];

//  epilog
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

#undef	p00_lo 
#undef	p00_hi 
#undef	p01_lo 
#undef	p01_hi 
#undef	p02_lo 
#undef	p02_hi 
#undef	p10_lo 
#undef	p10_hi 
#undef	p11_lo 
#undef	p11_hi 
#undef	p12_lo 
#undef	p12_hi 
#undef	p20_lo 
#undef	p20_hi	
#undef	p21_lo 
#undef	p21_hi 
#undef	p22_lo 
#undef	p22_hi 
#undef	tp0_lo 
#undef	tp0_hi 
#undef	tp3_lo 
#undef	tp3_hi 
#undef	tp2_lo 
#undef	tp2_hi 
#undef	tp1_lo 
#undef	tp1_hi

#undef	loop_i 
#undef	loop_o

#undef	Xlo 	
#undef	Xhi 	
#undef	Ylo 	
#undef	Yhi 	
#undef	Zlo 	
#undef	Zhi

#undef	bf0_lo 
#undef	bf0_hi 
#undef	bf1_lo 
#undef	bf1_hi 
#undef	bf2_lo 
#undef	bf2_hi 

    ret

    END