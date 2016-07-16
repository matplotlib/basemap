#!/bin/sh
area=10000000000000000000000000000000
gmt pscoast -M -Dc -R0/360/-90/90 -N1 -A${area} > countries_c.txt
gmt pscoast -M -Dc -R0/360/-90/90 -N2 -A${area} > states_c.txt
gmt pscoast -M -Dl -R0/360/-90/90 -N1 -A${area} > countries_l.txt
gmt pscoast -M -Dl -R0/360/-90/90 -N2 -A${area} > states_l.txt
gmt pscoast -M -Di -R0/360/-90/90 -N1 -A${area} > countries_i.txt
gmt pscoast -M -Di -R0/360/-90/90 -N2 -A${area} > states_i.txt
gmt pscoast -M -Dh -R0/360/-90/90 -N1 -A${area} > countries_h.txt
gmt pscoast -M -Dh -R0/360/-90/90 -N2 -A${area} > states_h.txt
gmt pscoast -M -Df -R0/360/-90/90 -N1 -A${area} > countries_f.txt
gmt pscoast -M -Df -R0/360/-90/90 -N2 -A${area} > states_f.txt
gmt pscoast -M -Dh -R0/360/-90/90 -N2 -A${area} > states_h.txt
gmt pscoast -M -Dc -R0/360/-90/90 -Ir -A${area} > rivers_c.txt
gmt pscoast -M -Dl -R0/360/-90/90 -Ir -A${area} > rivers_l.txt
gmt pscoast -M -Di -R0/360/-90/90 -Ir -A${area} > rivers_i.txt
gmt pscoast -M -Dh -R0/360/-90/90 -Ir -A${area} > rivers_h.txt
gmt pscoast -M -Df -R0/360/-90/90 -Ir -A${area} > rivers_f.txt
