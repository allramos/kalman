#include "contiki.h"
#include "contiki-conf.h"
#include "net/rime/rimeaddr.h"
#include "net/rime/collect.h"
#include "dev/cc2420.h"
#include "dev/cc2420_const.h"
#include <stdio.h> /* For printf() */
#include <stdlib.h>
#include "node-id.h"
#include <math.h>
#include <string.h>
#include <float.h>

#include "tinyekf_config.h"


		static int choldc1(double * a, double * p, int n) {
			int i,j,k;
			double sum;

			for (i = 0; i < n; i++) {
				for (j = i; j < n; j++) {
					sum = a[i*n+j];
					for (k = i - 1; k >= 0; k--) {
						sum -= a[i*n+k] * a[j*n+k];
					}
					if (i == j) {
						if (sum <= 0) {
                    return 1; /* error */
						}
                ///p[i] = sqrt(sum);
						p[i] = powf(sum, 0.5);
					}
					else {
						a[j*n+i] = sum / p[i];
					}
				}
			}

    return 0; /* success */
		}

		static int choldcsl(double * A, double * a, double * p, int n) 
		{
			int i,j,k; double sum;
			for (i = 0; i < n; i++) 
				for (j = 0; j < n; j++) 
					a[i*n+j] = A[i*n+j];
				if (choldc1(a, p, n)) return 1;
				for (i = 0; i < n; i++) {
					a[i*n+i] = 1 / p[i];
					for (j = i + 1; j < n; j++) {
						sum = 0;
						for (k = i; k < j; k++) {
							sum -= a[j*n+k] * a[k*n+i];
						}
						a[j*n+i] = sum / p[j];
					}
				}

    return 0; /* success */
			}


			static int cholsl(double * A, double * a, double * p, int n) 
			{
				int i,j,k;
				if (choldcsl(A,a,p,n)) return 1;
				for (i = 0; i < n; i++) {
					for (j = i + 1; j < n; j++) {
						a[i*n+j] = 0.0;
					}
				}
				for (i = 0; i < n; i++) {
					a[i*n+i] *= a[i*n+i];
					for (k = i + 1; k < n; k++) {
						a[i*n+i] += a[k*n+i] * a[k*n+i];
					}
					for (j = i + 1; j < n; j++) {
						for (k = j; k < n; k++) {
							a[i*n+j] += a[k*n+i] * a[k*n+j];
						}
					}
				}
				for (i = 0; i < n; i++) {
					for (j = 0; j < i; j++) {
						a[i*n+j] = a[j*n+i];
					}
				}

    return 0; /* success */
			}

			static void zeros(double * a, int m, int n)
			{
				int j;
				for (j=0; j<m*n; ++j)
					a[j] = 0;
			}

#ifdef DEBUG
			static void dump(double * a, int m, int n, const char * fmt)
			{
				int i,j;

				char f[100];
				sprintf(f, "%s ", fmt);
				for(i=0; i<m; ++i) {
					for(j=0; j<n; ++j)
						printf(f, a[i*n+j]);
					printf("\n");
				}
			}
#endif

/* C <- A * B */
			static void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols)
			{
				int i, j,l;

				for(i=0; i<arows; ++i)
					for(j=0; j<bcols; ++j) {
						c[i*bcols+j] = 0;
						for(l=0; l<acols; ++l)
							c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
					}
				}

				static void mulvec(double * a, double * x, double * y, int m, int n)
				{
					int i, j;

					for(i=0; i<m; ++i) {
						y[i] = 0;
						for(j=0; j<n; ++j)
							y[i] += x[j] * a[i*n+j];
					}
				}

				static void transpose(double * a, double * at, int m, int n)
				{
					int i,j;

					for(i=0; i<m; ++i)
						for(j=0; j<n; ++j) {
							at[j*m+i] = a[i*n+j];
						}
					}

/* A <- A + B */
					static void accum(double * a, double * b, int m, int n)
					{        
						int i,j;

						for(i=0; i<m; ++i)
							for(j=0; j<n; ++j)
								a[i*n+j] += b[i*n+j];
						}

/* C <- A + B */
						static void add(double * a, double * b, double * c, int n)
						{
							int j;

							for(j=0; j<n; ++j)
								c[j] = a[j] + b[j];
						}


/* C <- A - B */
						static void sub(double * a, double * b, double * c, int n)
						{
							int j;

							for(j=0; j<n; ++j)
								c[j] = a[j] - b[j];
						}

						static void negate(double * a, int m, int n)
						{        
							int i, j;

							for(i=0; i<m; ++i)
								for(j=0; j<n; ++j)
									a[i*n+j] = -a[i*n+j];
							}

							static void mat_addeye(double * a, int n)
							{
								int i;
								for (i=0; i<n; ++i)
									a[i*n+i] += 1;
							}

/* TinyEKF code ------------------------------------------------------------------- */

#include "tiny_ekf.h"

							typedef struct {

    double * x;    /* state vector */

    double * P;  /* prediction error covariance */
    double * Q;  /* process noise covariance */
    double * R;  /* measurement error covariance */

    double * G;  /* Kalman gain; a.k.a. K */

    double * F;  /* Jacobian of process model */
    double * H;  /* Jacobian of measurement model */

    double * Ht; /* transpose of measurement Jacobian */
    double * Ft; /* transpose of process Jacobian */
    double * Pp; /* P, post-prediction, pre-update */

    double * fx;  /* output of user defined f() state-transition function */
    double * hx;  /* output of user defined h() measurement function */

    /* temporary storage */
								double * tmp0;
								double * tmp1;
								double * tmp2;
								double * tmp3;
								double * tmp4;
								double * tmp5; 

							} ekf_t2;

							static void unpack(void * v, ekf_t2 * ekf, int n, int m)
							{
    /* skip over n, m in data structure */
								char * cptr = (char *)v;
								cptr += 2*sizeof(int);

								double * dptr = (double *)cptr;
								ekf->x = dptr;
								dptr += n;
								ekf->P = dptr;
								dptr += n*n;
								ekf->Q = dptr;
								dptr += n*n;
								ekf->R = dptr;
								dptr += m*m;
								ekf->G = dptr;
								dptr += n*m;
								ekf->F = dptr;
								dptr += n*n;
								ekf->H = dptr;
								dptr += m*n;
								ekf->Ht = dptr;
								dptr += n*m;
								ekf->Ft = dptr;
								dptr += n*n;
								ekf->Pp = dptr;
								dptr += n*n;
								ekf->fx = dptr;
								dptr += n;
								ekf->hx = dptr;
								dptr += m;
								ekf->tmp0 = dptr;
								dptr += n*n;
								ekf->tmp1 = dptr;
								dptr += n*m;
								ekf->tmp2 = dptr;
								dptr += m*n;
								ekf->tmp3 = dptr;
								dptr += m*m;
								ekf->tmp4 = dptr;
								dptr += m*m;
								ekf->tmp5 = dptr;
							}

							void ekf_init(void * v, int n, int m)
							{
    /* retrieve n, m and set them in incoming data structure */
								int * ptr = (int *)v;
								*ptr = n;
								ptr++;
								*ptr = m;

    /* unpack rest of incoming structure for initlization */
								ekf_t2 ekf;
								unpack(v, &ekf, n, m);

    /* zero-out matrices */
								zeros(ekf.P, n, n);
								zeros(ekf.Q, n, n);
								zeros(ekf.R, m, m);
								zeros(ekf.G, n, m);
								zeros(ekf.F, n, n);
								zeros(ekf.H, m, n);
}
double pPoints[2];
int ekf_step(void * v, double * z)
							{        
    /* unpack incoming structure */

								int * ptr = (int *)v;
								int n = *ptr;
								ptr++;
								int m = *ptr;

								ekf_t2 ekf;
								unpack(v, &ekf, n, m); 

    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
								mulmat(ekf.F, ekf.P, ekf.tmp0, n, n, n);
								transpose(ekf.F, ekf.Ft, n, n);
								mulmat(ekf.tmp0, ekf.Ft, ekf.Pp, n, n, n);
								accum(ekf.Pp, ekf.Q, n, n);
								




    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
								transpose(ekf.H, ekf.Ht, m, n);
								mulmat(ekf.Pp, ekf.Ht, ekf.tmp1, n, n, m);
								mulmat(ekf.H, ekf.Pp, ekf.tmp2, m, n, n);
								mulmat(ekf.tmp2, ekf.Ht, ekf.tmp3, m, n, m);
								accum(ekf.tmp3, ekf.R, m, m);
								if (cholsl(ekf.tmp3, ekf.tmp4, ekf.tmp5, m)) return 1;
								mulmat(ekf.tmp1, ekf.tmp4, ekf.G, n, m, m);

    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
								sub(z, ekf.hx, ekf.tmp5, m);
								mulvec(ekf.G, ekf.tmp5, ekf.tmp2, n, m);
								add(ekf.fx, ekf.tmp2, ekf.x, n);

    /* P_k = (I - G_k H_k) P_k */
								mulmat(ekf.G, ekf.H, ekf.tmp0, n, m, n);
								negate(ekf.tmp0, n, n);
								mat_addeye(ekf.tmp0, n);
								mulmat(ekf.tmp0, ekf.Pp, ekf.P, n, n, n);



    /* success */
								return 0;
							}



















//Trilat ==========================================
							int nos[3] = {-1,-1,-1};
							int sinais[3] = {0,0,0};

							static int TAM = 3;
							static int LIMIAR = -90;
//EKF ==========================================
							static const double T = 1;
							static const double veloc = 0.3;

							static int estados = 4;
							static int medicoes = 2;
//==========================================

							static struct collect_conn tc;
							static int cont = 1;

							ekf_t eKalman;

							PROCESS(example_collect_process, "RSS Measurement");
							AUTOSTART_PROCESSES(&example_collect_process);

							static double rssiToMeters(signed char r){
    signed char t; //rssi à um metro
    double A, B; //constantes
    double y; //saída

    t = -11;
    A = -0.587278856870205;
    B = 5.55352952962953;

    y = A * r - B;

    return y;
}

static int getSubString(char *source, char *target,int from, int to)
{
	int length=0;
	int i=0,j=0;
	
	//get length
	while(source[i++]!='\0')
		length++;
	
	if(from<0 || from>length){
		printf("Invalid \'from\' index\n");
		return 1;
	}
	if(to>length){
		printf("Invalid \'to\' index\n");
		return 1;
	}	
	
	for(i=from,j=0;i<=to;i++,j++){
		target[j]=source[i];
	}
	
	//assign NULL at the end of string
	target[j]='\0'; 
	
	return 0;	
}

typedef struct point{
	float x,y;
}Point;

float norm (Point p) // get the norm of a vector
{
	return powf(powf(p.x,2)+powf(p.y,2),.5);
}

Point trilateration(Point Point1, Point Point2, Point Point3, double r1, double r2, double r3) {
	Point position;

	double p2p1Distance = powf(powf(Point2.x-Point1.x,2) + powf(Point2.y - Point1.y,2),0.5);
	Point ex = {(Point2.x-Point1.x)/p2p1Distance, (Point2.y-Point1.y)/p2p1Distance};
	Point aux = {Point3.x-Point1.x,Point3.y-Point1.y};

	double i = ex.x * aux.x + ex.y * aux.y;

	Point aux2 = { Point3.x-Point1.x-i*ex.x, Point3.y-Point1.y-i*ex.y};
	Point ey = { aux2.x / norm (aux2), aux2.y / norm (aux2) };

	double j = ey.x * aux.x + ey.y * aux.y;

	double x = (powf(r1,2) - powf(r2,2) + powf(p2p1Distance,2))/ (2 * p2p1Distance);
	double y = (powf(r1,2) - powf(r3,2) + powf(i,2) + powf(j,2))/(2*j) - i*x/j;

	double finalX = Point1.x+ x*ex.x + y*ey.x;
	double finalY = Point1.y+ x*ex.y + y*ey.y;

	position.x = finalX;
	position.y = finalY;

	return position;
}

static int conectado(int origem, int *v, int t){
	int i;
	for(i = 0; i < t; i++){
		if(origem == v[i])
			return 1;
	}
	return 0;
}

static int conectados(int *v, int t){
	int cont = 0, i = 0;
	for(i = 0; i < t; i++){
		if(v[i] > 0)
			cont++;
	}
	return cont;
}

static int conexao(int origem, int sinal, Point pos, int *n, int *s, Point *p, int limiar, int t, uint8_t hops){
	int i;
	for(i = 0; i < t; i++){
		if(!conectado(origem, n, t) && hops == 1){
			if(n[i] == -1){ //novo nó. Entra com posição vazia
				n[i] = origem;
				s[i] = sinal;
				p[i].x = pos.x;
				p[i].y = pos.y;
			}					
		}else{
			if(origem == n[i])
				s[i] = sinal;
		}

		if(s[i] <= limiar){ //desconectar nós com limiar excedido
			n[i] = -1;
			s[i] = 0;
			p[i].x = 0;
			p[i].y = 0;
		}
	}
	return 0;
}
//=============================================================================

static void init(ekf_t * ekf, int xInicio, int yInicio) {
	double processNoise = 0.001; //0.1
	double measurementNoise = 0;
	double stateNoise = 0.1; //0.001
	int i;


	ekf->Q[0][0] = processNoise;
	ekf->Q[1][1] = processNoise;
	ekf->Q[2][2] = 0;
	ekf->Q[3][3] = 0;

	for (i = 0; i < estados; ++i)
		ekf->P[i][i] = stateNoise;

	for (i = 0; i < medicoes; ++i)
		ekf->R[i][i] = measurementNoise;

	//printf("Matriz Q = %d \n", (int) ekf->Q[0][0]);	

	// position
	ekf->x[0] = xInicio; //initial x
	ekf->x[1] = yInicio; //initial y

	// velocity
	ekf->x[2] = 0.1; //0
	ekf->x[3] = 0.1; //0

}

//static void model(ekf_t * ekf, double SV[2][1]) {
static void model(ekf_t * ekf) {
	int i, j;
	int T = 1;


	ekf->fx[0] = ekf->x[0];// + T * ekf->x[2];
	ekf->fx[1] = ekf->x[1];// + T * ekf->x[3];
	ekf->fx[2] = ekf->x[2];
	ekf->fx[3] = ekf->x[3];

	for (j = 0; j < estados; ++j)
		ekf->F[j][j] = 1;

	ekf->F[0][2] = 0.1; 
	ekf->F[1][3] = 0.1;
	

	ekf->hx[0] = ekf->x[0];// * 0.9;
	ekf->hx[1] = ekf->x[1];// * 0.9;

	
	ekf->H[0][0] = 1;
	// ekf->H[0][2] = 1; //Opção 1
	ekf->H[1][1] = 1;
	// ekf->H[1][3] = 1; //Opção 1
	ekf->H[2][2] = 1; //Opção 2
	ekf->H[3][3] = 1; //Opção 2
	
}


//=============================================================================
Point posicoes[3];
static void recv(const rimeaddr_t *originator, uint8_t seqno, uint8_t hops)
{
	static signed char rss;
	static signed char rss_val;
	static signed char rss_offset;  

	rss_val = cc2420_last_rssi;
	rss_offset=-45;
	rss=rss_val + rss_offset;

//Tratamento da msg contendo as coordenadas dos SNs
	char * msg = packetbuf_dataptr();

	char x[5], y[5];

	if(getSubString(msg,x,0,2) == 1)
		getSubString(msg,x,0,3);

	if(getSubString(msg,y,3, strlen(msg) - 1) == 1)
		getSubString(msg,y,4, strlen(msg) - 1);

	int origem = originator->u8[0];	
	int pos_x = atoi(x);
	int pos_y = atoi(y);
	Point posicao;
	posicao.x = pos_x;
	posicao.y = pos_y;
//====================================================
//INICIO - Trilateração
	
	Point finalPos;

	conexao(origem, rss, posicao, nos, sinais, posicoes, LIMIAR, TAM, hops);

	//printf("No: %d RSSI: %d x: %d y: %d     ", origem, rss, pos_x, pos_y);
	printf("[ %d | %d | %d ]  ", nos[0], nos[1], nos[2]);
	if(conectados(nos, TAM) <3)
		printf("\n");

	if(conectados(nos, TAM) == 3){
		Point p1 = {posicoes[0].x,posicoes[0].y};
		Point p2 = {posicoes[1].x,posicoes[1].y};
		Point p3 = {posicoes[2].x,posicoes[2].y};

		double d1 = rssiToMeters(sinais[0]);
		double d2 = rssiToMeters(sinais[1]);
		double d3 = rssiToMeters(sinais[2]);

		finalPos = trilateration(p1,p2,p3,d1,d2,d3);
		

		//INICIO - Extended Kalman Filter =========================================

		double res[4];
		res[0] = eKalman.x[0];
		res[1] = eKalman.x[1];
		res[2] = eKalman.x[2];
		res[3] = eKalman.x[3];
		int i, j, k;
		for(i = 0; i < 1; i ++){
			
				for (k=0; k<4; k++){ // 2
					res[k] += eKalman.fx[k] * eKalman.x[k];
			 	} // 2
			
		}




		double z[medicoes];
		z[0] = (int) finalPos.x;
		z[1] = (int) finalPos.y;

		// z[0] = (int) node_loc_x;
		// z[1] = (int) node_loc_y;

		model(&eKalman);

		ekf_step(&eKalman, z);

		double xFinal = eKalman.x[0];
		double yFinal = eKalman.x[1];
		// ===========================
		
		
		// ===========================
		printf("Real  ->  x: %d   y: %d\n", (int) node_loc_x, (int) node_loc_y);

		printf("               Trilat->  x: %d   y: %d\n", (int) finalPos.x, (int) finalPos.y);

		printf("               Filtro->  x: %f   y: %f\n", xFinal, yFinal);		

		// printf("               Predic->  x: %d   y: %d\n\n", (int) res[0], (int) res[1]);

		//FIM    - Extended Kalman Filter ==============================

	}//end if
	
//FIM    - Trilateração ==============================

}



static const struct collect_callbacks callbacks = { recv };

PROCESS_THREAD(example_collect_process, ev, data)
{
	static struct etimer periodic;
	static struct etimer et;
	static int contador = 0;
	static int continua = 1;

	PROCESS_BEGIN();

	collect_open(&tc, 130, COLLECT_ROUTER, &callbacks);

	if(rimeaddr_node_addr.u8[0] == 1 && rimeaddr_node_addr.u8[1] == 0) 
	{
		printf("I am sink\n");
		collect_set_sink(&tc, 1);


		ekf_init(&eKalman, 4, 2);

			// Do local initialization
		init(&eKalman, 0, 0);
	}

  /* Allow some time for the network to settle. */
  etimer_set(&et, 12 * CLOCK_SECOND);//120
  //etimer_set(&et, 10 * CLOCK_SECOND);
  PROCESS_WAIT_UNTIL(etimer_expired(&et));

  //extendedKalmanFilter(); //================================== INICIO EKF

  //=======================================================================

  while(continua == 1) 
  {
  //	if(node_loc_x > 50)
    //    	continua = 0;
    /* Send a packet every 1 seconds. */
  	if(etimer_expired(&periodic)) 
  	{
  		etimer_set(&periodic, CLOCK_SECOND * 1 );
      //etimer_set(&et, random_rand() % (CLOCK_SECOND * 1)); //Enviar em tempo aleatório entre 0 e 1s
      etimer_set(&et, (CLOCK_SECOND * 1)); //Enviar a cada 1s
  }

  PROCESS_WAIT_EVENT();


  if(etimer_expired(&et)) 
  {
  	static rimeaddr_t oldparent;
  	const rimeaddr_t *parent;
  	if(rimeaddr_node_addr.u8[0] != 1 )
  	{

  		packetbuf_clear();
        //packetbuf_set_datalen(sprintf(packetbuf_dataptr(),"%s", "Fight On") + 1);
  		packetbuf_set_datalen(sprintf(packetbuf_dataptr(), "%d %d", node_loc_x, node_loc_y) + 1);
        //packetbuf_set_datalen(sprintf(packetbuf_dataptr(),"%s", "Oi") + 1);
  		contador++;

  		collect_send(&tc, 15);

  		parent = collect_parent(&tc);
  		if(!rimeaddr_cmp(parent, &oldparent)) 
  		{
  			if(!rimeaddr_cmp(&oldparent, &rimeaddr_null))
  			{
  				printf("#L %d 0\n", oldparent.u8[0]);
  			}
  			if(!rimeaddr_cmp(parent, &rimeaddr_null)) 
  			{
  				printf("#L %d 1\n", parent->u8[0]);
  			}

  			rimeaddr_copy(&oldparent, parent);
  		}        
  	}
  }

  } //end of while
  //printf("FINAL");
  PROCESS_END();
} //end of process thread




