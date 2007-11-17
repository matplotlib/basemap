#!/bin/sh
area=10000000000000000000000000000000
pscoast -M -Dc -R0/360/-90/90 -N1 -JQ180/4.5i -W -A${area} > countries_c.txt
pscoast -M -Dc -R0/360/-90/90 -N2 -JQ180/4.5i -W -A${area} > states_c.txt
pscoast -M -Dl -R0/360/-90/90 -N1 -JQ180/4.5i -W -A${area} > countries_l.txt
pscoast -M -Dl -R0/360/-90/90 -N2 -JQ180/4.5i -W -A${area} > states_l.txt
pscoast -M -Di -R0/360/-90/90 -N1 -JQ180/4.5i -W -A${area} > countries_i.txt
pscoast -M -Di -R0/360/-90/90 -N2 -JQ180/4.5i -W -A${area} > states_i.txt
pscoast -M -Dh -R0/360/-90/90 -N1 -JQ180/4.5i -W -A${area} > countries_h.txt
pscoast -M -Dh -R0/360/-90/90 -N2 -JQ180/4.5i -W -A${area} > states_h.txt
pscoast -M -Df -R0/360/-90/90 -N1 -JQ180/4.5i -W -A${area} > countries_f.txt
pscoast -M -Df -R0/360/-90/90 -N2 -JQ180/4.5i -W -A${area} > states_f.txt
pscoast -M -Dh -R0/360/-90/90 -N2 -JQ180/4.5i -W -A${area} > states_h.txt
pscoast -M -Dc -R0/360/-90/90 -Ir -JQ180/4.5i -W -A${area} > rivers_c.txt
pscoast -M -Dl -R0/360/-90/90 -Ir -JQ180/4.5i -W -A${area} > rivers_l.txt
pscoast -M -Di -R0/360/-90/90 -Ir -JQ180/4.5i -W -A${area} > rivers_i.txt
pscoast -M -Dh -R0/360/-90/90 -Ir -JQ180/4.5i -W -A${area} > rivers_h.txt
pscoast -M -Df -R0/360/-90/90 -Ir -JQ180/4.5i -W -A${area} > rivers_f.txt
