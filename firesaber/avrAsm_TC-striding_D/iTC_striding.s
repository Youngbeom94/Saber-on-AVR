            
            NAME    iTC_striding_asm     ;   module name
            PUBLIC  iTC_striding_asm     ;   make the main label vissible outside this module
            RSEG    CODE                 ;   place program in 'CODE' segment            
      
mc_mul16x16_lo  MACRO m_rlo, m_rhi, m_s1lo, m_s1hi, m_s2lo, m_s2hi

	mul     m_s1lo, m_s2lo
    movw    m_rlo,  r0
    mul     m_s1lo, m_s2hi
    add     m_rhi,  r0
    mul     m_s1hi, m_s2lo
    add     m_rhi,  r0
    
    ENDM
    
mc_mul16x8_lo  MACRO m_rlo, m_rhi, m_s1lo, m_s1hi, m_s2lo

	mul     m_s1lo, m_s2lo
    movw    m_rlo,  r0 
    mul     m_s1hi, m_s2lo
    add     m_rhi,  r0
    
    ENDM

mc_mul16x8_lo_acc  MACRO m_rlo, m_rhi, m_s1lo, m_s1hi, m_s2lo

	mul     m_s1lo, m_s2lo
    add     m_rlo,  r0 
    adc     m_rhi,  r1
    mul     m_s1hi, m_s2lo
    add     m_rhi,  r0
    
    ENDM

mc_mul16x8_lo_scc  MACRO m_rlo, m_rhi, m_s1lo, m_s1hi, m_s2lo

	mul     m_s1lo, m_s2lo
    sub     m_rlo,  r0 
    sbc     m_rhi,  r1
    mul     m_s1hi, m_s2lo
    sub     m_rhi,  r0
    
    ENDM   
   
   
iTC_striding_asm:

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
                                
#define s0_lo       r2	        
#define s0_hi       r3	        
#define s1_lo       r4	        
#define s1_hi       r5	        
#define s2_lo       r6	        
#define s2_hi       r7	        
#define s3_lo       r8	        
#define s3_hi       r9	        
#define s4_lo       r10	        
#define s4_hi       r11	        
#define s5_lo       r12	        
#define s5_hi       r13	        
#define s6_lo       r14	        
#define s6_hi       r15	        
#define t0_lo       r16	        
#define t0_hi       r17	        
#define t1_lo       r18	        
#define t1_hi       r19
#define t2_lo       r20	        
#define t2_hi       r21	
#define t3_lo       r26	        
#define t3_hi       r27	

#define zero		r22	
#define tmp         r23		
                                                                                              
#define ofs         r24	        
#define loop_o      r25

#define Xlo         r26	        
#define Xhi         r27	        
#define Ylo         r28	        
#define Yhi         r29
#define Zlo         r30	        
#define Zhi         r31


initialize:
    
    movw    Zlo,    r16 ; Y = Z = src1_ex
    movw    Ylo,    r16 

    clr     zero
    ldi     ofs,    128    
    ldi     loop_o, 64    
                
LOOP:

    // Load src
    ldd     s0_lo,  Z+0
    ldd     s0_hi,  Z+1
    add     Zlo,    ofs
    adc     Zhi,    zero
    
    ldd     s1_lo,  Z+0
    ldd     s1_hi,  Z+1
    add     Zlo,    ofs
    adc     Zhi,    zero
    
    ldd     s2_lo,  Z+0
    ldd     s2_hi,  Z+1
    add     Zlo,    ofs
    adc     Zhi,    zero
    
    ldd     s3_lo,  Z+0
    ldd     s3_hi,  Z+1
    add     Zlo,    ofs
    adc     Zhi,    zero
    
    ldd     s4_lo,  Z+0
    ldd     s4_hi,  Z+1
    add     Zlo,    ofs
    adc     Zhi,    zero
    
    ldd     s5_lo,  Z+0
    ldd     s5_hi,  Z+1
    add     Zlo,    ofs
    adc     Zhi,    zero
    
    ldd     s6_lo,  Z+0
    ldd     s6_hi,  Z+1
    add     Zlo,    ofs
    adc     Zhi,    zero
    
    // {-360, -120,  -40,    5,     3,    8, -360}
    movw    t0_lo,  s0_lo
    add     t0_lo,  s6_lo
    adc     t0_hi,  s6_hi
    movw    t1_lo,  t0_lo
    lsl     t0_lo
    rol     t0_hi
    add     t0_lo,  t1_lo
    adc     t0_hi,  t1_hi   ; t0 = (res[0][i] + res[6][i]) * 3

    add     t0_lo,  s1_lo
    adc     t0_hi,  s1_hi   ; t0 = (res[0][i] + res[6][i]) * 3 + res[1][i]

    movw    t1_lo,  t0_lo
    lsl     t0_lo
    rol     t0_hi
    add     t0_lo,  t1_lo
    adc     t0_hi,  t1_hi   ; t0 = ((res[0][i] + res[6][i]) * 3 + res[1][i] ) * 3

    add     t0_lo,  s2_lo
    adc     t0_hi,  s2_hi   ; t0 = ( (res[0][i] + res[6][i]) * 3 + res[1][i]) * 3 + res[2][i];
        
    movw    t2_lo,  t0_lo
    lsl     t2_lo
    rol     t2_hi
    lsl     t2_lo
    rol     t2_hi
    add     t2_lo,  t0_lo
    adc     t2_hi,  t0_hi   ; t2 = t0 * 5
    
    movw    t0_lo,  s5_lo  
    sub     t0_lo,  t2_lo
    sbc     t0_hi,  t2_hi   ; t0 = res[5][i] - t0 * 5;
    
    lsl     t0_lo
    rol     t0_hi
    lsl     t0_lo
    rol     t0_hi
    lsl     t0_lo
    rol     t0_hi           ; t0 *= 8

    movw    t2_lo,  s3_lo
    lsl     t2_lo
    rol     t2_hi
    lsl     t2_lo
    rol     t2_hi
    add     t2_lo,  s3_lo
    adc     t2_hi,  s3_hi
    add     t0_lo,  t2_lo
    adc     t0_hi,  t2_hi   ; t0 += res[3][i] * 5;

    movw    t2_lo,  s4_lo    
    lsl     t2_lo
    rol     t2_hi
    add     t2_lo,  s4_lo
    adc     t2_hi,  s4_hi
    add     t0_lo,  t2_lo
    adc     t0_hi,  t2_hi   ; t0 += res[4][i] * 3;    
    
    ldi     t1_lo,  0xA5
    ldi     t1_hi,  0x4F        
    mc_mul16x16_lo  t2_lo, t2_hi, t0_lo, t0_hi, t1_lo, t1_hi    ; t2 = TC4_buff[1] = t0 * iTC4_modified_scale[1]; // 0x4FA5    
    
    lsr     t2_hi
    ror     t2_lo
    lsr     t2_hi
    ror     t2_lo   ; t2 = TC4_buff[1] >> 2;

    push    t2_lo                                                                                           ; TC4_buff[1]_lo -> stack 1
    push    t2_hi                                                                                           ; TC4_buff[1]_hi -> stack 2

    // { -30,   16,   16,   -1,    -1,    0,   96}
    movw    t0_lo,  s1_lo
    add     t0_lo,  s2_lo
    adc     t0_hi,  s2_hi   ; t0 = res[1][i] + res[2][i];
    
    movw    t1_lo,  s3_lo
    add     t1_lo,  s4_lo
    adc     t1_hi,  s4_hi   ; t1 = res[3][i] + res[4][i];

    lsl     t0_lo
    rol     t0_hi
    lsl     t0_lo
    rol     t0_hi
    lsl     t0_lo
    rol     t0_hi
    lsl     t0_lo
    rol     t0_hi
    sub     t0_lo,  t1_lo
    sbc     t0_hi,  t1_hi   ; t0 = 16 * t0 - t1;

    ldi     t1_lo,  30
    mc_mul16x8_lo  t2_lo, t2_hi, s0_lo, s0_hi, t1_lo
    sub     t0_lo,  t2_lo
    sbc     t0_hi,  t2_hi   ; t0 = t0 - 30 * res[0][i];

    ldi     t1_lo,  96
    mc_mul16x8_lo_acc  t0_lo, t0_hi, s6_lo, s6_hi, t1_lo    ; t0 += 96 * res[6][i];
    
    ldi     t1_lo,  0xAB
    ldi     t1_hi,  0xAA        
    mc_mul16x16_lo  t2_lo, t2_hi, t0_lo, t0_hi, t1_lo, t1_hi    ; t2 = TC4_buff[1] = t0 * iTC4_modified_scale[1]; // 0x4FA5 

    lsr     t2_hi
    ror     t2_lo
    lsr     t2_hi
    ror     t2_lo
    lsr     t2_hi
    ror     t2_lo   ; t2 = TC4_buff[2] >> 3;
    
    push    t2_lo                                                                                           ; TC4_buff[2]_lo -> stack 3
    push    t2_hi                                                                                           ; TC4_buff[2]_hi -> stack 4
    
    // {  45,   27,   -7,   -1,     0,   -1,   45},
    movw    t0_lo,  s0_lo
    add     t0_lo,  s6_lo
    adc     t0_hi,  s6_hi   ; t0 = res[0][i] + res[6][i];

    movw    t1_lo,  s3_lo
    add     t1_lo,  s5_lo
    adc     t1_hi,  s5_hi   ; t1 = res[3][i] + res[5][i];

    ldi     t3_lo,  45
    mc_mul16x8_lo  t2_lo, t2_hi, t0_lo, t0_hi, t3_lo
    sub     t2_lo,  t1_lo
    sbc     t2_hi,  t1_hi
    movw    t0_lo,  t2_lo   ; t0 = t2 = 45 * t0 - t1;

    ldi     t1_lo,  27
    mc_mul16x8_lo_acc  t0_lo, t0_hi, s1_lo, s1_hi, t1_lo    ; t0 += 27 * res[1][i];

    ldi     t3_lo,  7
    mc_mul16x8_lo_scc  t0_lo, t0_hi, s2_lo, s2_hi, t3_lo    ; t0 = t0 - 7 * res[2][i];            
    
    ldi     t1_lo,  0x39
    ldi     t1_hi,  0x8E        
    mc_mul16x16_lo  t2_lo, t2_hi, t0_lo, t0_hi, t1_lo, t1_hi    ; t2 = TC4_buff[3] = t0 * iTC4_modified_scale[3]; // 0x8E39
    
    lsr     t2_hi
    ror     t2_lo   ; t2 = TC4_buff[3] >> 1;
    
    push    t2_lo                                                                                           ; TC4_buff[3]_lo -> stack 5
    push    t2_hi                                                                                           ; TC4_buff[3]_hi -> stack 6

    // {   6,   -4,   -4,    1,     1,    0, -120}
    movw    t0_lo,  s1_lo
    add     t0_lo,  s2_lo
    adc     t0_hi,  s2_hi   ; t0 = res[1][i] + res[2][i];

    movw    t1_lo,  s3_lo
    add     t1_lo,  s4_lo
    adc     t1_hi,  s4_hi   ; t1 = res[3][i] + res[4][i];

    lsl     t0_lo
    rol     t0_hi
    lsl     t0_lo
    rol     t0_hi
    sub     t1_lo,  t0_lo
    sbc     t1_hi,  t0_hi
    movw    t0_lo,  t1_lo   ; t0 = t1 - 4 * t0;

    ldi     t1_lo,  6
    mc_mul16x8_lo_acc  t0_lo, t0_hi, s0_lo, s0_hi, t1_lo    ; t0 += 27 * res[1][i];

    ldi     t3_lo,  120
    mc_mul16x8_lo_scc  t0_lo, t0_hi, s6_lo, s6_hi, t3_lo    ; t0 = t0 - 7 * res[2][i];    

    ldi     t1_lo,  0xAB
    ldi     t1_hi,  0xAA        
    mc_mul16x16_lo  t2_lo, t2_hi, t0_lo, t0_hi, t1_lo, t1_hi    ; t2 = TC4_buff[4] = t0 * iTC4_modified_scale[4];  // 0xAAAB
    
    lsr     t2_hi
    ror     t2_lo
    lsr     t2_hi
    ror     t2_lo
    lsr     t2_hi
    ror     t2_lo   ; t2 = TC4_buff[4] >> 3;
    
    push    t2_lo                                                                                           ; TC4_buff[4]_lo -> stack 7
    push    t2_hi                                                                                           ; TC4_buff[4]_hi -> stack 8

    // { -90,  -60,   20,    5,    -3,    2,  -90},
    movw    t0_lo,  s0_lo
    add     t0_lo,  s6_lo
    adc     t0_hi,  s6_hi   ; t0 = res[0][i] + res[6][i];    

    ldi     t1_lo,  0xA6
    ldi     t1_hi,  0xFF 
    mc_mul16x16_lo  t2_lo, t2_hi, t0_lo, t0_hi, t1_lo, t1_hi    ; t0 *= -90;
    movw    t0_lo,  t2_lo

    ldi     t3_lo,  60
    mc_mul16x8_lo_scc  t0_lo, t0_hi, s1_lo, s1_hi, t3_lo    ; t0 = t0 - 60 * res[1][i];

    ldi     t1_lo,  20
    mc_mul16x8_lo_acc  t0_lo, t0_hi, s2_lo, s2_hi, t1_lo    ; t0 += 20 * res[2][i];

    ldi     t1_lo,  5
    mc_mul16x8_lo_acc  t0_lo, t0_hi, s3_lo, s3_hi, t1_lo    ; t0 += 5 * res[3][i];

    ldi     t3_lo,  3
    mc_mul16x8_lo_scc  t0_lo, t0_hi, s4_lo, s4_hi, t3_lo    ; t0 = t0 - 3 * res[4][i];

    lsl     s5_lo
    rol     s5_hi
    add     t0_lo, s5_lo
    adc     t0_hi, s5_hi
    
    ldi     t1_lo,  0xA5
    ldi     t1_hi,  0x4F        
    mc_mul16x16_lo  t2_lo, t2_hi, t0_lo, t0_hi, t1_lo, t1_hi    ; t2 = TC4_buff[5] = t0 * iTC4_modified_scale[5]; // 0x4FA5
            
    lsr     t2_hi
    ror     t2_lo
    lsr     t2_hi
    ror     t2_lo   ; t2 = TC4_buff[4] >> 3;
    
    push    t2_lo                                                                                           ; TC4_buff[5]_lo -> stack 9
    push    t2_hi                                                                                           ; TC4_buff[5]_hi -> stack 10 
    
    // Multiply by powers of two.    
    pop     s5_hi                                                                                           ; stack -> TC4_buff[5]_hi 9
    pop     s5_lo                                                                                           ; stack -> TC4_buff[5]_lo 8
    
    pop     s4_hi                                                                                           ; stack -> TC4_buff[4]_hi 7
    pop     s4_lo                                                                                           ; stack -> TC4_buff[4]_lo 6
    
    pop     s3_hi                                                                                           ; stack -> TC4_buff[3]_hi 5
    pop     s3_lo                                                                                           ; stack -> TC4_buff[3]_lo 4
    
    pop     s2_hi                                                                                           ; stack -> TC4_buff[2]_hi 3
    pop     s2_lo                                                                                           ; stack -> TC4_buff[2]_lo 2
        
    pop     s1_hi                                                                                           ; stack -> TC4_buff[1]_hi 1
    pop     s1_lo                                                                                           ; stack -> TC4_buff[1]_lo 0

    movw    Zlo,    Ylo
    add     Zlo,    ofs
    adc     Zhi,    zero
    std     Z+0,    s1_lo
    std     Z+1,    s1_hi   ;   res[1][i] = TC4_buff[1] >> 2;
    
    add     Zlo,    ofs
    adc     Zhi,    zero
    std     Z+0,    s2_lo
    std     Z+1,    s2_hi   ;   res[2][i] = TC4_buff[2] >> 3;

    add     Zlo,    ofs
    adc     Zhi,    zero
    std     Z+0,    s3_lo
    std     Z+1,    s3_hi   ;   res[3][i] = TC4_buff[3] >> 1;

    add     Zlo,    ofs
    adc     Zhi,    zero
    std     Z+0,    s4_lo
    std     Z+1,    s4_hi   ;   res[3][i] = TC4_buff[4] >> 3;

    add     Zlo,    ofs
    adc     Zhi,    zero
    std     Z+0,    s5_lo
    std     Z+1,    s5_hi   ;   res[3][i] = TC4_buff[4] >> 3;
                   
    adiw    Ylo, 2       ; offset control
    movw    Zlo,    Ylo    
    dec     loop_o
    breq    LOOP_END
    rjmp    LOOP
    
LOOP_END:

#undef s0_lo             
#undef s0_hi             
#undef s1_lo             
#undef s1_hi             
#undef s2_lo             
#undef s2_hi             
#undef s3_lo             
#undef s3_hi             
#undef s4_lo     	        
#undef s4_hi     	        
#undef s5_lo     	        
#undef s5_hi     	        
#undef s6_lo     	        
#undef s6_hi     	        
#undef t0_lo     	        
#undef t0_hi     	        
#undef t1_lo     	        
#undef t1_hi     
#undef t2_lo     	        
#undef t2_hi     	
#undef t3_lo     	        
#undef t3_hi     	
               
#undef zero		
#undef tmp       		
                                                                                     
#undef ofs       	        
#undef loop_o    
               
#undef Xlo       	        
#undef Xhi       	        
#undef Ylo       	        
#undef Yhi       
#undef Zlo       	        
#undef Zhi       

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

    ret

    END