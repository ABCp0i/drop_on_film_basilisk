#rm _post.c post
#serial
#qcc -Wall -g -C -grid=octree  -std=c99 post.c -o post -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

#openmp 
qcc -Wall -fopenmp -std=c99 -O2 post.c -o post_op -I$HOME -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU  -lm -lOSMesa 

#mpi
#qcc -source -grid=octree -D_dropstats=0 -D_MPI=1 post.c
#mpicc -Wall -std=c99 -O2 _post.c -o post_mpi -I$BASILISK  -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
