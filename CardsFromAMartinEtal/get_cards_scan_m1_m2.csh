#! /bin/tcsh

set m1 = 150
while ($m1 <= 1000)
  set m2 = 150
  while ($m2 <= 1000)
    sed "s/.*M_1.*/        1    ${m1}    # M_1/" suspect2_lha.in.in > suspect2_lha.in.temp
    sed "s/.*M_2.*/        2    ${m2}    # M_2/" suspect2_lha.in.temp > suspect2_lha.in
    ./run
    mv susyhit_slha.out cards/susyhit_slha_mu_m150_M1_${m1}_M2_${m2}.out
    @ m2 = $m2 + 50
  end
  @ m1 = $m1 + 50
end
