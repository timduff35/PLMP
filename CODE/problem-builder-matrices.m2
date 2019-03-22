load "problem-builder-core.m2"
-- set up relevant matrices only

pointMatrices = apply(lineIncidences, LINES->commonPointMatrix(C,L/(l->l^LINES)))
lineMatrices = apply(nLines, i->seeLineMatrix(C,L/(l->l^{i})))
end
