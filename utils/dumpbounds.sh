#!/bin/sh
set -e

# This should be set to wherever you extract the GSHHG data to.
gshhs_dir=./gshhg/
area=10000000000000000000000000000000

gmt gshhg $gshhs_dir/gshhs_c.b > gshhs_c.txt
gmt gshhg $gshhs_dir/gshhs_l.b > gshhs_l.txt
gmt gshhg $gshhs_dir/gshhs_i.b > gshhs_i.txt
gmt gshhg $gshhs_dir/gshhs_h.b > gshhs_h.txt
gmt gshhg $gshhs_dir/gshhs_f.b > gshhs_f.txt

gmt pscoast -M -Rg -Dc -N1 -A${area} > countries_c.txt
gmt pscoast -M -Rg -Dc -N2 -A${area} > states_c.txt
gmt pscoast -M -Rg -Dl -N1 -A${area} > countries_l.txt
gmt pscoast -M -Rg -Dl -N2 -A${area} > states_l.txt
gmt pscoast -M -Rg -Di -N1 -A${area} > countries_i.txt
gmt pscoast -M -Rg -Di -N2 -A${area} > states_i.txt
gmt pscoast -M -Rg -Dh -N1 -A${area} > countries_h.txt
gmt pscoast -M -Rg -Dh -N2 -A${area} > states_h.txt
gmt pscoast -M -Rg -Df -N1 -A${area} > countries_f.txt
gmt pscoast -M -Rg -Df -N2 -A${area} > states_f.txt
gmt pscoast -M -Rg -Dc -Ir -A${area} > rivers_c.txt
gmt pscoast -M -Rg -Dl -Ir -A${area} > rivers_l.txt
gmt pscoast -M -Rg -Di -Ir -A${area} > rivers_i.txt
gmt pscoast -M -Rg -Dh -Ir -A${area} > rivers_h.txt
gmt pscoast -M -Rg -Df -Ir -A${area} > rivers_f.txt
