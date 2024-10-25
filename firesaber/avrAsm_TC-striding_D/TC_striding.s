            
            NAME    TC_striding_asm     ;   module name
            PUBLIC  TC_striding_asm     ;   make the main label vissible outside this module
            RSEG    CODE                ;   place program in 'CODE' segment            
   
TC_striding_asm:

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

#define zero  	r1

#define s0_lo 	r2  
#define s0_hi 	r3  
#define s1_lo 	r4  
#define s1_hi 	r5  
#define s2_lo 	r6  
#define s2_hi 	r7  
#define s3_lo 	r8  
#define s3_hi 	r9  
#define t0_lo 	r10 
#define t0_hi 	r11 
#define t1_lo 	r12 
#define t1_hi 	r13 
#define t2_lo 	r14 
#define t2_hi 	r15 
#define t3_lo 	r16 
#define t3_hi 	r17 
#define tp0_lo	r18 
#define tp0_hi	r19 

#define ofs     r22                
#define loop_o	r23
#define adr_lo  r24
#define adr_hi  r25

#define Xlo		r26
#define Xhi		r27
#define Ylo		r28
#define Yhi		r29

initialize:
    
    movw    Ylo,    r16 ; Y = src1_ex
    movw    Xlo,    r20 ; X = src

    clr     zero
    ldi     ofs,    128    
    ldi     loop_o, 64
    movw    adr_lo, Ylo
                
LOOP:
    
    ld      s0_lo,  X+
    ld      s0_hi,  X+
    
    ld      s1_lo,  X+
    ld      s1_hi,  X+
    
    ld      s2_lo,  X+
    ld      s2_hi,  X+
    
    ld      s3_lo,  X+
    ld      s3_hi,  X+
    
    std     Y+0,    s0_lo   ; src1_extended[0][i] = src1[i * 4 + 0];
    std     Y+1,    s0_hi
    
    movw    t0_lo,  s0_lo
    add     t0_lo,  s2_lo
    adc     t0_hi,  s2_hi   ; t0 = src1[i * 4 + 0] + src1[i * 4 + 2];
    
    movw    t1_lo,  s1_lo
    add     t1_lo,  s3_lo
    adc     t1_hi,  s3_hi   ; t1 = src1[i * 4 + 1] + src1[i * 4 + 3];
    
    movw    tp0_lo, t0_lo   ; tp0 = t0
    add     tp0_lo, t1_lo
    adc     tp0_hi, t1_hi   ; tp0 = t0 + t1

    sub     t0_lo,  t1_lo
    sbc     t0_hi,  t1_hi   ; t0 = t0 - t1
    
    add     Ylo,    ofs
    adc     Yhi,    zero
    std     Y+0,    tp0_lo   
    std     Y+1,    tp0_hi   ; src1_extended[1][i] = t0 + t1;

    add     Ylo,    ofs
    adc     Yhi,    zero
    std     Y+0,    t0_lo   
    std     Y+1,    t0_hi   ; src1_extended[2][i] = t0 - t1;


    movw    t2_lo,  s2_lo
    lsl     t2_lo
    rol     t2_hi
    lsl     t2_lo
    rol     t2_hi
    add     t2_lo,  s0_lo
    adc     t2_hi,  s0_hi   ; t2 = src1[i * 4 + 0] + 4 * src1[i * 4 + 2];
    
    movw    t3_lo,  s3_lo
    lsl     t3_lo
    rol     t3_hi
    lsl     t3_lo
    rol     t3_hi
    add     t3_lo,  s1_lo
    adc     t3_hi,  s1_hi   ; t3 = src1[i * 4 + 1] + 4 * src1[i * 4 + 3];
    lsl     t3_lo
    rol     t3_hi           ; t3 = 2 * t3;

    movw    tp0_lo, t2_lo   ; tp0 = t2
    add     tp0_lo, t3_lo
    adc     tp0_hi, t3_hi   ; tp0 = t2 + 2 * t3

    sub     t2_lo, t3_lo
    sbc     t2_hi, t3_hi    ; t2 - 2 * t3;

    add     Ylo,    ofs
    adc     Yhi,    zero
    std     Y+0,    tp0_lo   
    std     Y+1,    tp0_hi  ; src1_extended[3][i] = t2 + 2 * t3;

    add     Ylo,    ofs
    adc     Yhi,    zero
    std     Y+0,    t2_lo   
    std     Y+1,    t2_hi   ; src1_extended[4][i] = t2 - 2 * t3;
        
    lsl     s0_lo
    rol     s0_hi
    add     s0_lo,  s1_lo
    adc     s0_hi,  s1_hi
    lsl     s0_lo
    rol     s0_hi
    add     s0_lo,  s2_lo
    adc     s0_hi,  s2_hi
    lsl     s0_lo
    rol     s0_hi
    add     s0_lo,  s3_lo
    adc     s0_hi,  s3_hi
    
    add     Ylo,    ofs
    adc     Yhi,    zero
    std     Y+0,    s0_lo   
    std     Y+1,    s0_hi   ; src1_extended[5][i] = ((src1[i * 4 + 0] * 2 + src1[i * 4 + 1]) * 2 + src1[i * 4 + 2] * 2 + src1[i * 4 + 3];

    add     Ylo,    ofs
    adc     Yhi,    zero
    std     Y+0,    s3_lo   
    std     Y+1,    s3_hi   ; src1_extended[6][i] = src1[i * 4 + 3];
        
    adiw    adr_lo, 2       ; offset control
    movw    Ylo,    adr_lo    
    dec     loop_o
    breq    LOOP_END
    rjmp    LOOP
    
LOOP_END:


#undef zero  
               
#undef s0_lo  
#undef s0_hi  
#undef s1_lo  
#undef s1_hi  
#undef s2_lo  
#undef s2_hi  
#undef s3_lo  
#undef s3_hi  
#undef t0_lo  
#undef t0_hi  
#undef t1_lo  
#undef t1_hi  
#undef t2_lo  
#undef t2_hi  
#undef t3_lo  
#undef t3_hi  
#undef tp0_lo 
#undef tp0_hi 
               
#undef ofs                    
#undef loop_o
#undef adr_lo  
#undef adr_hi  
               
#undef Xlo	
#undef Xhi	
#undef Ylo	
#undef Yhi	


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