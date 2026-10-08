#! /bin/sh
if [ ! -e done ]; then
    cd ../Cards/
    ln -sf ../GrnFnc_blvblv/*_card.dat .
    cd ../SubProcesses/
    cd ./P1_gg_mu+vmmu-vmxbbx/
    ln -sf ../../GrnFnc_blvblv/*.tbl .
    ln -sf ../../GrnFnc_blvblv/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvblv/matrix_1_gg_mu+vmmu-vmxbbx.f ./matrix.f
    ln -sf ../../GrnFnc_blvblv/leshouche_sng.inc leshouche.inc
    cd ../P8_gg_mu+vmmu-vmxbbx/
    ln -sf ../../GrnFnc_blvblv/*.tbl .
    ln -sf ../../GrnFnc_blvblv/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvblv/matrix_8_gg_mu+vmmu-vmxbbx.f ./matrix.f
    cd ../P8_uux_mu+vmmu-vmxbbx/
    ln -sf ../../GrnFnc_blvblv/*.tbl .
    ln -sf ../../GrnFnc_blvblv/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvblv/matrix_8_uux_mu+vmmu-vmxbbx.f ./matrix.f
    cd ../P8_uxu_mu+vmmu-vmxbbx/
    ln -sf ../../GrnFnc_blvblv/*.tbl .
    ln -sf ../../GrnFnc_blvblv/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvblv/matrix_8_uxu_mu+vmmu-vmxbbx.f ./matrix.f
    cd ../P8_ddx_mu+vmmu-vmxbbx/
    ln -sf ../../GrnFnc_blvblv/*.tbl .
    ln -sf ../../GrnFnc_blvblv/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvblv/matrix_8_ddx_mu+vmmu-vmxbbx.f ./matrix.f
    cd ../P8_dxd_mu+vmmu-vmxbbx/
    ln -sf ../../GrnFnc_blvblv/*.tbl .
    ln -sf ../../GrnFnc_blvblv/{thr,kfac,coupling,kin_functions}.inc .
    mv matrix.f matrix_org.f
    ln -sf ../../GrnFnc_blvblv/matrix_8_dxd_mu+vmmu-vmxbbx.f ./matrix.f
    cd ../../GrnFnc_blvblv/
    touch done
else
    echo "done"
fi
