computeDimDegreeF4 = method(Options=>{RandomSeed=>0,Field=>ZZ/nextPrime 10000,FileName=>null,Saturate=>false})
computeDimDegreeF4 (ZZ,ZZ,List) := o -> (l,g,i) -> (
    D = (l,g,i);
    setRandomSeed o.RandomSeed;
    FF = o.Field;
    f = if o.FileName === null then stdio 
    else if fileExists o.FileName then return
    else openOut o.FileName;
    f << "D = " << D << endl; 
    setupTime := first elapsedTiming load "problem-builder.m2";
    f << "setup time = " << setupTime << " sec" << endl;
    f << "F4 time = " << first elapsedTiming (gbJz := groebnerBasis(Jz, Strategy=>"F4")) << " sec"<< endl;
    inJ := if o.Saturate then (
	J := ideal apply(flatten entries selectInSubring(1,gbJz), g->sub(g,R));
	scan(lineIncidences, LINES -> J = saturate(J,commonPointSat(C,L/(l->l^LINES))));
	scan(nLines, i -> J = saturate(J,seeLineSat(C,L/(l->l^{i}))));
	ideal leadTerm J
	) else ideal apply(flatten entries selectInSubring(1,gbJz), g->sub(leadMonomial g,R));
    if inJ == 0 then f << "zero ideal!!!" << endl else (  
    	f << "dim = " << dim inJ << endl;
    	f << "deg = " << degree inJ << endl;
    	);
    f << close
    ) 

isMinimal = method()
isMinimal (ZZ,ZZ,List) := (l,g,i) -> (
    D = (l,g,i);
    isParametric = true;
    elapsedTime load "problem-builder-matrices.m2";
    xy := fabricatedVectorFLPQ D;
    assert all(pointMatrices,M->rank sub(M, xy)<4);
    assert all(lineMatrices,M->rank sub(M, xy)<3);
    matrices := pointMatrices | lineMatrices;
    gm := goodMinors(pointMatrices,lineMatrices,xy);
    << "D = "; print D;
    << "goodMINORS = "; print gm;
--    << toString(#gm) << " is the numbr of good minors" << endl;
    return (#gm == (m-1)*6-1);
    )
