# Serial
#qcc -Wall -grid=quadtree -std=c99 -DFILTERED=1 -O2 impact.c -o impact -I$HOME -L$HOME/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

# openMP
qcc -Wall -fopenmp -grid=quadtree -std=c99 -DFILTERED=1 -O2 impact.c -o impact -I$HOME -L$HOME/gl -lm
