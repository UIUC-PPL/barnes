./charmrun +p4 ./barnes -killat=4 -dtime=0.00001 -ppc=100 -p=20 -b=8 -in=1k.dat -balancePeriod=2 +balancer Orb3dLB_notopo +noAnytimeMigration +LBPeriod 0.01 +LBCommOff > out.lb
#for i in *.dot; do dot -Tps $i -o lb.`basename $i .dot`.ps; done
#rm *.dot
./charmrun +p4 ./barnes -killat=4 -dtime=0.00001 -ppc=100 -p=20 -b=8 -in=1k.dat -balancePeriod=4 +balancer Orb3dLB_notopo +noAnytimeMigration +LBPeriod 0.01 +LBCommOff > out.nolb
#for i in *.dot; do dot -Tps $i -o nolb.`basename $i .dot`.ps; done
#rm *.dot
echo "sorting output"
grep "^(.*,.*) " out.lb | sort -n > sorted.out.lb
grep "^(.*,.*) " out.nolb | sort -n > sorted.out.nolb

echo "output accelerations "
grep final sorted.out.lb > final.lb
grep final sorted.out.nolb > final.nolb
