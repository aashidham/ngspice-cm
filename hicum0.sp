
Hicum0 test .option
          temp=25.0 reltol=1e-5 bstol=1e-15 
          .model mybjtmodel npn
          npn=1 level=7 qN1 c b 0 0 mybjtmodel rcc cc c 100.0 rbb
          bb b 1k vcc cc 0 1.0 vbb bb 0 0.0 ac=0.0 .dc vbb 0 1.0
          0.02 .print dc v(c) i(vbb) i(vcc) .end
