#!/usr/bin/env python

import sys
import string
file=sys.argv[1]
enabled = ["ISS","PO101","SO50","AO91","LILACSAT2","AO27"]
whitelist = string.ascii_letters + string.digits 
with open(file,"r") as rd, open("generated/sat_tles.c","w") as wd:
        lines = rd.readlines()
        for i in range(0,len(lines),3):
            item="".join(x for x in lines[i].strip() if x in whitelist)
            print(item)
            line1 = lines[i+1].strip()
            line2 = lines[i+2].strip()
            if item in ["ISS","PO101","SO50","AO91","LILACSAT2","AO27"]:
                j = enabled.index(item)
                s="""
                parse_elements( 
                \"%s\", 
                \"%s\", 
                &(satellites[%d].tle));
                strncpy(satellites[%d].name, \"%s\", 15);
                satellites[%d].name[15] = 0;
                num_satellites += 1;
                """%(line1, line2, j, j, item, j)
                wd.write(s);
                i+=1


# PO101 schedule: https://nitter.net/Diwata2PH
