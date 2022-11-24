set datafile separator ','

set xlabel 'Time'
d(y) = ($0 == 0) ? (y1 = y, 1/0) : (y2 = y1, y1 = y, y1-y2)
d2(x,y) = ($0 == 0) ? (x1 = x, y1 = y, 1/0) : (x2 = x1, x1 = x, y2 = y1, y1 = y, (y1-y2)/(x1-x2))
set ylabel "degrees"
#set y2tics
#set ytics nomirror
#set y2label "meters"
#set y2range [0:]
#set y2range [0:*] reverse

set style line 100 lw 0.5 lt rgb "#000000" 
set style line 101 lw 0.1 lt rgb "#ff0000" 
set style line 102 lw 0.5 lt rgb "#50a050" 
set style line 103 lw 1 lt rgb "#b8ee30" 
set style line 104 lw 1 lt rgb "#00ffff" 
set style line 105 lw 1 lt rgb "#ffff00" 
set style line 106 lw 2 lt rgb "#ff00ff" 

#plot "out.csv" using 1:2 with lines ls 101 title "az", '' using 1:3 with lines ls 102 title "elev", '' using 1:4 with lines ls 103 title "dist" axis x1y2 
#plot 'out.csv' using 1:3 with lines ls 102 title "elev", '' using 1:4 with lines ls 103 title "dist" axis x1y2, 0 ls 104
#plot 'out.csv' using 1:3 with lines ls 102 title "elev", 0 ls 104
#plot 'out.csv' using 1:3 with lines ls 102 title "elev",  0 ls 104, '' using 1:(d($3)) title '1-variable derivative' ls 101
#plot 'out.csv' using 1:3 with lines ls 102 title "elev",  0 ls 104, '' using 1:(d2($1,$3)) title '1-variable derivative' ls 101

#plot 'ISS.csv' using 1:3 with lines ls 102 title "elev", 0 ls 100, \
#'ISS_test.csv' using 1:3 with dots ls 106 title 'elevtest',\
#'ISS_test.csv' using 1:3 with lines ls 101 title 'elevtestl'
#'ISS_test.csv' using 1:(d2($1,$3)) title 'elevtest d' ls 100

#plot 'ISS.csv' using 1:3 with lines ls 101 title "ISS", 0 ls 100, \
#'PO101.csv' using 1:3 with lines ls 102 title "PO101", \
#'SO50.csv' using 1:3 with lines ls 103 title "SO50", \
#'AO91.csv' using 1:3 with lines ls 104 title "AO91", \
#'LILACSAT2.csv' using 1:3 with lines ls 105 title "LILACSAT2", \
#'AO27.csv' using 1:3 with lines ls 106 title "AO27"


plot '25544_groundtruth.csv' using 1:3 with lines ls 102 title "elev", \
0 ls 100, \
'25544_mmdev.csv' using 1:3 with dots ls 106 title 'elevtest',\
'25544_mmdev.csv' using 1:3 with lines ls 101 title 'elevtestl'
