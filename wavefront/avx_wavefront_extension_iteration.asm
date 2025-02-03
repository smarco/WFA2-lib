section .data
vecShuffle db 60,61,62,63,56,57,58,59,52,53,54,55,48,49,50,51,44,45,46,47,40,41,42,43,36,37,38,39,32,33,34,35,28,29,30,31,24,25,26,27,20,21,22,23,16,17,18,19,12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3

section .text
global avx_wavefront_extension_iteration
global load_avx2_sequence
default rel

avx_wavefront_extension_iteration:
    xor rax, rax

    vpxord zmm6, zmm6, zmm6
    vpandnd zmm5, zmm6, zmm6
    vmovdqu8 zmm10, [vecShuffle]

    vmovdqu32 zmm0, [rdi]
    vmovdqu32 zmm1, zmm0
    vmovdqu32 zmm2, [rsi]
    vmovdqu32 zmm3, [rdx]
    vpsubd zmm4, zmm0, zmm2

    vmovdqu32 [rcx], zmm4
    vpaddd zmm4, zmm2, zmm3

    vmovdqu32 [rsi], zmm4

    vpcmpgtd k1, zmm0, zmm5

    vpgatherdd zmm7{k1}, [r8+zmm4*1]
    vpgatherdd zmm8{k1}, [r9+zmm1*1]

    vpcmpeqd k2, zmm7, zmm8

    vpxord zmm9, zmm7, zmm8
    vpshufb zmm9, zmm9, zmm10

    vplzcntd zmm11{k1}{z}, zmm9
    vpsrld zmm12, zmm11, 3
    vpaddd zmm0{k1}{z}, zmm0, zmm12

    vmovdqu32 [rdi], zmm0

    ret


load_avx2_sequence:
    ; assignment
    vmovdqu ymm0, [rdi]
    vmovdqu ymm1, [rsi]

    ; compare sequences
    vpcmpeqb ymm2, ymm0, ymm1

    ; return to third parameter
    vmovdqu [rdx], ymm2
    vmovdqu [rsi], ymm0

    vpcmpeqb ymm2, ymm0, ymm1



    ; first pass
    vextracti128 xmm3, ymm2, 0
    vextracti128 xmm4, ymm2, 1
    vpaddb xmm3, xmm3, xmm4

    ; second pass
    vpextrq r8, xmm3, 0
    vpextrq r9, xmm3, 1
    movq xmm4, r8
    movq xmm5, r9
    vpaddb xmm3, xmm4, xmm5

    pxor xmm4, xmm4
    pxor xmm5, xmm5

    ; third pass
    vpextrd rcx, xmm3, 0
    vpextrd rdx, xmm3, 1
    movq xmm4, rcx
    movq xmm5, rdx
    vpaddb xmm3, xmm4, xmm5

    pxor xmm4, xmm4
    pxor xmm5, xmm5
    xor ecx, ecx
    xor edx, edx

    ; fourth pass
    vpextrw rcx, xmm3, 0
    vpextrw rdx, xmm3, 1
    movq xmm4, rcx
    movq xmm5, rdx
    vpaddb xmm3, xmm4, xmm5

    pxor xmm4, xmm4
    pxor xmm5, xmm5
    xor rcx, rcx
    xor rdx, rdx

    ; fifth pass
    vpextrb rcx, xmm3, 0
    vpextrb rdx, xmm3, 1
    movq xmm4, rcx
    movq xmm5, rdx
    vpaddb xmm3, xmm4, xmm5

    mov eax, 0xFF
    vmovd ebx, xmm3
    sub eax, ebx
    add eax, 0x01

    ret