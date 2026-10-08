#! /bin/sh
if [ ! -e done ]; then
    cd ../Cards/
    ln -sf ../GrnFnc_bwbw/*_card.dat .
    cd ../SubProcesses/
    cd ./P1_gg_bw+bxw-/
    ln -sf ../../GrnFnc_bwbw/*.tbl .
    ln -sf ../../GrnFnc_bwbw/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_bwbw/matrix_1_gg_bw+bxw-.f ./matrix.f 
    ln -sf ../../GrnFnc_bwbw/leshouche_sng.inc leshouche.inc
    cd ../P8_gg_bw+bxw-/
    ln -sf ../../GrnFnc_bwbw/*.tbl .
    ln -sf ../../GrnFnc_bwbw/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_bwbw/matrix_8_gg_bw+bxw-.f ./matrix.f
    cd ../P8_uux_bw+bxw-/
    ln -sf ../../GrnFnc_bwbw/*.tbl .
    ln -sf ../../GrnFnc_bwbw/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_bwbw/matrix_8_uux_bw+bxw-.f ./matrix.f
    cd ../P8_uxu_bw+bxw-/
    ln -sf ../../GrnFnc_bwbw/*.tbl .
    ln -sf ../../GrnFnc_bwbw/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_bwbw/matrix_8_uxu_bw+bxw-.f ./matrix.f
    cd ../P8_ddx_bw+bxw-/
    ln -sf ../../GrnFnc_bwbw/*.tbl .
    ln -sf ../../GrnFnc_bwbw/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_bwbw/matrix_8_ddx_bw+bxw-.f ./matrix.f
    cd ../P8_dxd_bw+bxw-/
    ln -sf ../../GrnFnc_bwbw/*.tbl .
    ln -sf ../../GrnFnc_bwbw/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_bwbw/matrix_8_dxd_bw+bxw-.f ./matrix.f
    cd ../../GrnFnc_bwbw/
    touch done
else
    echo "done"
fi
