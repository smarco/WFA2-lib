section .data
vecShuffle db 60,61,62,63,56,57,58,59,52,53,54,55,48,49,50,51,44,45,46,47,40,41,42,43,36,37,38,39,32,33,34,35,28,29,30,31,24,25,26,27,20,21,22,23,16,17,18,19,12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3

section .text
global avx_wavefront_extension_iteration
global load_avx2_sequence
default rel

avx_wavefront_extension_iteration:
    xor rax, rax

    vpxord zmm6, zmm6, zmm6
    vpternlogd zmm5, zmm6, zmm6, 0xFF
    vmovdqu8 zmm10, [vecShuffle]

    vmovdqu32 zmm0, [rdi]
    vmovdqu32 zmm1, zmm0
    vmovdqu32 zmm2, [rsi]
    vpsubd zmm4, zmm0, zmm2

    vpcmpgtd k1, zmm0, zmm6
    kmovd k3, k1

    vpxord zmm7, zmm7, zmm7
    vpgatherdd zmm7{k3}, [rcx+zmm4*1]
    kmovd k3, k1

    vpxord zmm8, zmm8, zmm8
    vpgatherdd zmm8{k3}, [r8+zmm1*1]

    vpcmpeqd k2{k1}, zmm7, zmm8

    kmovd [rdx], k2

    vpxord zmm9, zmm7, zmm8
    vpshufb zmm9, zmm9, zmm10

    vmovdqu32 [r9], zmm9

    vplzcntd zmm11{k1}{z}, zmm9
    
    vpsrld zmm12, zmm11, 3

    vpaddd zmm12, zmm0, zmm12
    vpxord zmm0, zmm0, zmm0
    vmovdqa32 zmm0{k1}{z}, zmm12

    sub rsp, 16
    mov rax, rsp
    add rsp, 16

    mov dword [rax], 0x10
    vpbroadcastd zmm3, [rax]
    vpaddd zmm4, zmm2, zmm3
    vmovdqu32 [rsi], zmm4

    vmovdqu32 [rdi], zmm12

    ret