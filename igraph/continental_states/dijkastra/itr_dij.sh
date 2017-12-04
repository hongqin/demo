 gmn3 states.config states.dic _gmn3.config 1  
 network_2_pairs.0.0.pl _gmn3.config _pairs.tab    
 ccom02 _pairs.tab _ccom02.rpt states.dic _largest  
 dijkastra_by_categories.0.0.pl -i _largest -r _dij.largest.rpt -q 1   
 cat _dij.largest.rpt >> dij.gmn3.rpt

