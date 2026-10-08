#! /bin/sh
if [ ! -e done ]; then
    cd ../Cards/
    ln -sf ../GrnFnc_blvbjj/*_card.dat .
    cd ../SubProcesses/
    cd ./P1_gg_mu+vmduxbbx/
    ln -sf ../../GrnFnc_blvbjj/*.tbl .
    ln -sf ../../GrnFnc_blvbjj/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvbjj/matrix_1_gg_mu+vmduxbbx.f ./matrix.f
    ln -sf ../../GrnFnc_blvbjj/leshouche_sng.inc leshouche.inc
    cd ../P8_gg_mu+vmduxbbx/
    ln -sf ../../GrnFnc_blvbjj/*.tbl .
    ln -sf ../../GrnFnc_blvbjj/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvbjj/matrix_8_gg_mu+vmduxbbx.f ./matrix.f
    cd ../P8_uux_mu+vmduxbbx/
    ln -sf ../../GrnFnc_blvbjj/*.tbl .
    ln -sf ../../GrnFnc_blvbjj/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvbjj/matrix_8_uux_mu+vmduxbbx.f ./matrix.f
    cd ../P8_uxu_mu+vmduxbbx/
    ln -sf ../../GrnFnc_blvbjj/*.tbl .
    ln -sf ../../GrnFnc_blvbjj/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvbjj/matrix_8_uxu_mu+vmduxbbx.f ./matrix.f
    cd ../P8_ddx_mu+vmduxbbx/
    ln -sf ../../GrnFnc_blvbjj/*.tbl .
    ln -sf ../../GrnFnc_blvbjj/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvbjj/matrix_8_ddx_mu+vmduxbbx.f ./matrix.f
    cd ../P8_dxd_mu+vmduxbbx/
    ln -sf ../../GrnFnc_blvbjj/*.tbl .
    ln -sf ../../GrnFnc_blvbjj/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvbjj/matrix_8_dxd_mu+vmduxbbx.f ./matrix.f
    cd ../../GrnFnc_blvbjj/
    touch done
else
    echo "done"
fi
