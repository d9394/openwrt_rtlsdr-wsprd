/*
 * FreeBSD License
 * Copyright (c) 2016, Guenael
 * All rights reserved.
 *
 * This file is based on rtl-sdr project, contribution :
 *   Copyright (C) 2012 by Steve Markgraf <steve@steve-m.de>
 *   Copyright (C) 2012 by Hoernchen <la@tfc-server.de>
 *   Copyright (C) 2012 by Kyle Keen <keenerd@gmail.com>
 *   Copyright (C) 2013 by Elias Oenal <EliasOenal@gmail.com>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <curl/curl.h>
#include <pthread.h>
#include <rtl-sdr.h>

#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <fftw3.h>
#include <sched.h>

#include "rtlsdr_wsprd.h"
#include "fano.h"
#include "nhash.h"
#include "wsprd_utils.h"
#include "wsprsim_utils.h"
#include "metric_tables.h"

/* TODO
 - multi device selection option
 - multispot report in one post
 - verbose option
*/

#define SIGNAL_LENGHT       116
#define SIGNAL_SAMPLE_RATE  375
#define SAMPLING_RATE       2400000
#define FS4_RATE            SAMPLING_RATE / 4                   // = 600 kHz
#define DOWNSAMPLING        SAMPLING_RATE / SIGNAL_SAMPLE_RATE  // = 6400
#define DEFAULT_BUF_LENGTH  (4 * 16384)                         // = 65536

#define DF       375.0/256.0
#define DT       1.0/375.0
#define TWOPIDT  2 * M_PI * DT
#define NIQ      45000
#define NBITS    81
#define NSYM     162
#define NSPERSYM 256
#define NFILT    256
#define NSIG     NSYM * NSPERSYM

/* Possible PATIENCE options: FFTW_ESTIMATE, FFTW_ESTIMATE_PATIENT, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE */
#define PATIENCE FFTW_ESTIMATE

uint8_t pr3[NSYM]= {
    1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0
};

fftwf_plan PLAN1,
           PLAN2,
           PLAN3;
int32_t printdata=0;

/* Global declaration for these structs */
struct receiver_state   rx_state;
struct receiver_options rx_options;
struct decoder_options  dec_options;
struct decoder_results  dec_results[50];
static rtlsdr_dev_t *rtl_device = NULL;


/* Thread stuff for separate decoding */
struct decoder_state {
    pthread_t        thread;
    pthread_attr_t   tattr;

    pthread_rwlock_t rw;
    pthread_cond_t   ready_cond;
    pthread_mutex_t  ready_mutex;
};
struct decoder_state dec;


/* Thread stuff for separate RX (blocking function) */
struct dongle_state {
    pthread_t        thread;
};
struct dongle_state dongle;

//***************************************************************************
void sync_and_demodulate(float *id, float *qd, long np,
                         uint8_t *symbols, float *f1, float fstep,
                         int32_t *shift1, int32_t lagmin, int32_t lagmax, int32_t lagstep,
                         float *drift1, int32_t symfac, float *sync, int32_t mode) {
    /***********************************************************************
     * mode = 0: no frequency or drift search. find best time lag.          *
     *        1: no time lag or drift search. find best frequency.          *
     *        2: no frequency or time lag search. calculate soft-decision   *
     *           symbols using passed frequency and shift.                  *
     ************************************************************************/

    float fbest=0.0;
    float f0=0.0,fp,ss;
    int32_t lag;
    static float fplast=-10000.0;
    float i0[NSYM],q0[NSYM],
          i1[NSYM],q1[NSYM],
          i2[NSYM],q2[NSYM],
          i3[NSYM],q3[NSYM];
    float p0,p1,p2,p3,cmet,totp,syncmax,fac;
    float c0[NSPERSYM],s0[NSPERSYM],
          c1[NSPERSYM],s1[NSPERSYM],
          c2[NSPERSYM],s2[NSPERSYM],
          c3[NSPERSYM],s3[NSPERSYM];
    float dphi0, cdphi0, sdphi0,
          dphi1, cdphi1, sdphi1,
          dphi2, cdphi2, sdphi2,
          dphi3, cdphi3, sdphi3;
    float fsum=0.0, f2sum=0.0, fsymb[NSYM];
    int32_t best_shift = 0;
    int32_t ifmin=0, ifmax=0;

    syncmax=-1e30;
    if( mode == 0 ) {
        ifmin=0;
        ifmax=0;
        fstep=0.0;
        f0=*f1;
    }
    if( mode == 1 ) {
        lagmin=*shift1;
        lagmax=*shift1;
        ifmin=-5;
        ifmax=5;
        f0=*f1;
    }
    if( mode == 2 ) {
        lagmin=*shift1;
        lagmax=*shift1;
        ifmin=0;
        ifmax=0;
        f0=*f1;
    }

    printf("Start sync and demodulate\n");

    for(int32_t ifreq=ifmin; ifreq<=ifmax; ifreq++) {
        f0=*f1+ifreq*fstep;
        for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) {
            ss=0.0;
            totp=0.0;
            for (int32_t i=0; i<NSYM; i++) {
                fp = f0 + ((float)*drift1/2.0)*((float)i-81.0)/81.0;
                if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
                    dphi0=TWOPIDT*(fp-1.5*DF);
                    cdphi0=cosf(dphi0);
                    sdphi0=sinf(dphi0);

                    dphi1=TWOPIDT*(fp-0.5*DF);
                    cdphi1=cosf(dphi1);
                    sdphi1=sinf(dphi1);

                    dphi2=TWOPIDT*(fp+0.5*DF);
                    cdphi2=cosf(dphi2);
                    sdphi2=sinf(dphi2);

                    dphi3=TWOPIDT*(fp+1.5*DF);
                    cdphi3=cosf(dphi3);
                    sdphi3=sinf(dphi3);

                    c0[0]=1;
                    s0[0]=0;
                    c1[0]=1;
                    s1[0]=0;
                    c2[0]=1;
                    s2[0]=0;
                    c3[0]=1;
                    s3[0]=0;

                    for (int32_t j=1; j<NSPERSYM; j++) {
                        c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                        s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                        c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                        s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                        c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                        s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                        c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                        s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
                    }
                    fplast = fp;
                }

                i0[i]=0.0;
                q0[i]=0.0;
                i1[i]=0.0;
                q1[i]=0.0;
                i2[i]=0.0;
                q2[i]=0.0;
                i3[i]=0.0;
                q3[i]=0.0;

                for (int32_t j=0; j<NSPERSYM; j++) {
                    int32_t k=lag+i*NSPERSYM+j;
                    if( (k>0) & (k<np) ) {
                        i0[i]=i0[i] + id[k]*c0[j] + qd[k]*s0[j];
                        q0[i]=q0[i] - id[k]*s0[j] + qd[k]*c0[j];
                        i1[i]=i1[i] + id[k]*c1[j] + qd[k]*s1[j];
                        q1[i]=q1[i] - id[k]*s1[j] + qd[k]*c1[j];
                        i2[i]=i2[i] + id[k]*c2[j] + qd[k]*s2[j];
                        q2[i]=q2[i] - id[k]*s2[j] + qd[k]*c2[j];
                        i3[i]=i3[i] + id[k]*c3[j] + qd[k]*s3[j];
                        q3[i]=q3[i] - id[k]*s3[j] + qd[k]*c3[j];
                    }
                }
                p0=i0[i]*i0[i] + q0[i]*q0[i];
                p1=i1[i]*i1[i] + q1[i]*q1[i];
                p2=i2[i]*i2[i] + q2[i]*q2[i];
                p3=i3[i]*i3[i] + q3[i]*q3[i];

                p0=sqrtf(p0);
                p1=sqrtf(p1);
                p2=sqrtf(p2);
                p3=sqrtf(p3);

                totp=totp+p0+p1+p2+p3;
                cmet=(p1+p3)-(p0+p2);
                ss=ss+cmet*(2*pr3[i]-1);
                if( mode == 2) {                 //Compute soft symbols
                    if(pr3[i]) {
                        fsymb[i]=p3-p1;
                    } else {
                        fsymb[i]=p2-p0;
                    }
                }
            }

            if( ss/totp > syncmax ) {          //Save best parameters
                syncmax=ss/totp;
                best_shift=lag;
                fbest=f0;
            }
        } // lag loop
    } //freq loop

    if( mode <=1 ) {                       //Send best params back to caller
        *sync=syncmax;
        *shift1=best_shift;
        *f1=fbest;
        return;
    }

    if( mode == 2 ) {
        *sync=syncmax;
        for (int32_t i=0; i<NSYM; i++) {              //Normalize the soft symbols
            fsum=fsum+fsymb[i]/NSYM;
            f2sum=f2sum+fsymb[i]*fsymb[i]/NSYM;
        }
        fac=sqrtf(f2sum-fsum*fsum);
        for (int32_t i=0; i<NSYM; i++) {
            fsymb[i]=symfac*fsymb[i]/fac;
            if( fsymb[i] > 127) fsymb[i]=127.0;
            if( fsymb[i] < -128 ) fsymb[i]=-128.0;
            symbols[i]=fsymb[i] + 128;
        }
        return;
    }
    return;
}


/***************************************************************************
 symbol-by-symbol signal subtraction
 ****************************************************************************/
void subtract_signal(float *id, float *qd, long np,
                     float f0, int32_t shift0, float drift0, uint8_t* channel_symbols) {

    float i0,q0;
    float c0[NSPERSYM],s0[NSPERSYM];
    float dphi, cdphi, sdphi;

    for (int32_t i=0; i<NSYM; i++) {
        float fp = f0 + ((float)drift0/2.0)*((float)i-81.0)/81.0;

        dphi=TWOPIDT*(fp+((float)channel_symbols[i]-1.5)*DF);
        cdphi=cosf(dphi);
        sdphi=sinf(dphi);

        c0[0]=1;
        s0[0]=0;

        for (int32_t j=1; j<NSPERSYM; j++) {
            c0[j]=c0[j-1]*cdphi - s0[j-1]*sdphi;
            s0[j]=c0[j-1]*sdphi + s0[j-1]*cdphi;
        }

        i0=0.0;
        q0=0.0;

        for (int32_t j=0; j<NSPERSYM; j++) {
            int32_t k=shift0+i*NSPERSYM+j;
            if( (k>0) & (k<np) ) {
                i0=i0 + id[k]*c0[j] + qd[k]*s0[j];
                q0=q0 - id[k]*s0[j] + qd[k]*c0[j];
            }
        }


        // subtract the signal here.

        i0=i0/NSPERSYM; //will be wrong for partial symbols at the edges...
        q0=q0/NSPERSYM;

        for (int32_t j=0; j<NSPERSYM; j++) {
            int32_t k=shift0+i*NSPERSYM+j;
            if( (k>0) & (k<np) ) {
                id[k]=id[k]- (i0*c0[j] - q0*s0[j]);
                qd[k]=qd[k]- (q0*c0[j] + i0*s0[j]);
            }
        }
    }
    return;
}


/******************************************************************************
 Fully coherent signal subtraction
 *******************************************************************************/
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int32_t shift0, float drift0, uint8_t* channel_symbols) {

    float phi=0, dphi, cs;
    float refi[NIQ]= {0}, refq[NIQ]= {0},
                                     ci[NIQ]= {0},   cq[NIQ]= {0},
                                             cfi[NIQ]= {0},  cfq[NIQ]= {0};

    /******************************************************************************
     Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
     Reference is:                       r(t) = exp( j*phi(t) )
     Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
     so c(t) has phase angle theta-phi
     Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
     *******************************************************************************/

    // create reference wspr signal vector, centered on f0.
    //
    for (int32_t i=0; i<NSYM; i++) {

        cs=(float)channel_symbols[i];

        dphi=TWOPIDT * ( f0 +
                         ((float)drift0/2.0)*((float)i-(float)NSYM/2.0)/((float)NSYM/2.0) +
                         (cs-1.5)*DF  );

        for (int32_t j=0; j<NSPERSYM; j++ ) {
            int32_t ii=NSPERSYM*i+j;
            refi[ii]=refi[ii]+cosf(phi); //cannot precompute sin/cos because dphi is changing
            refq[ii]=refq[ii]+sinf(phi);
            phi=phi+dphi;
        }
    }

    // s(t) * conjugate(r(t))
    // beginning of first symbol in reference signal is at i=0
    // beginning of first symbol in received data is at shift0.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (int32_t i=0; i<NSYM*NSPERSYM; i++) {
        int32_t k=shift0+i;
        if( (k>0) & (k<np) ) {
            ci[i+NFILT] = id[k]*refi[i] + qd[k]*refq[i];
            cq[i+NFILT] = qd[k]*refi[i] - id[k]*refq[i];
        }
    }

    //quick and dirty filter - may want to do better
    float w[NFILT]= {0}, norm=0, partialsum[NFILT]= {0};
    for (int32_t i=0; i<NFILT; i++) {
        w[i]=sinf(M_PI*(float)i/(float)(NFILT-1));
        norm=norm+w[i];
    }
    for (int32_t i=0; i<NFILT; i++) {
        w[i]=w[i]/norm;
    }
    for (int32_t i=1; i<NFILT; i++) {
        partialsum[i]=partialsum[i-1]+w[i];
    }

    // LPF
    for (int32_t i=NFILT/2; i<NIQ-NFILT/2; i++) {
        cfi[i]=0.0;
        cfq[i]=0.0;
        for (int32_t j=0; j<NFILT; j++) {
            cfi[i]=cfi[i]+w[j]*ci[i-NFILT/2+j];
            cfq[i]=cfq[i]+w[j]*cq[i-NFILT/2+j];
        }
    }

    // subtract c(t)*r(t) here
    // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq)+cq*refi)
    // beginning of first symbol in reference signal is at i=NFILT
    // beginning of first symbol in received data is at shift0.
    for (int32_t i=0; i<NSIG; i++) {
        if( i<NFILT/2 ) {        // take care of the end effect (LPF step response) here
            norm=partialsum[NFILT/2+i];
        } else if( i>(NSIG-1-NFILT/2) ) {
            norm=partialsum[NFILT/2+NSIG-1-i];
        } else {
            norm=1.0;
        }
        int32_t k=shift0+i;
        int32_t j=i+NFILT;
        if( (k>0) & (k<np) ) {
            id[k]=id[k] - (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
            qd[k]=qd[k] - (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
        }
    }
    return;
}

//***************************************************************************

int32_t wspr_decode(float *idat, float *qdat, uint32_t npoints,
                    struct decoder_options options, struct decoder_results *decodes,
                    int32_t *n_results) {

    printf("Start wspr decode\n");

    int32_t i,j,k;
    uint32_t metric, maxcycles, cycles, maxnp;
    uint8_t symbols[NBITS*2]= {0};
    uint8_t decdata[(NBITS+7)/8]= {0};
    int8_t message[]= {-9,13,-35,123,57,-39,64,0,0,0,0};

    char callsign[13]= {0};
    char call_loc_pow[23]= {0};
    char call[13]= {0};
    char loc[7]= {0};
    char pwr[3]= {0};

    int32_t delta,verbose=0;
    int32_t writenoise=0,wspr_type=2, ipass;
    int32_t shift1, lagmin, lagmax, lagstep, worth_a_try, not_decoded;
    float freq0[200],snr0[200],drift0[200],sync0[200];
    int32_t shift0[200];
    float dt_print;
    double freq_print;
    double dialfreq= (double)options.freq / 1e6; // check
    float dialfreq_error=0.0;
    float f1, fstep, sync1=0.0, drift1;
    int32_t noprint=0;
    int32_t uniques=0;
    float fmin=-110.0;
    float fmax=110.0;
    char hashtab[32768*13]= {0};
    int32_t nh;

    float allfreqs[100];
    char allcalls[100][13];
    
    memset(allfreqs,0,sizeof(float)*100);
    memset(allcalls,0,sizeof(char)*100*13);

    printf("memset ready\n");

    // Parameters used for performance-tuning:
    maxcycles=10000;                         //Fano timeout limit
    double minsync1=0.10;                    //First sync limit
    double minsync2=0.12;                    //Second sync limit
    int32_t iifac=3;                             //Step size in final DT peakup
    int32_t symfac=50;                           //Soft-symbol normalizing factor
    int32_t maxdrift=4;                          //Maximum (+/-) drift
    double minrms=52.0 * (symfac/64.0);      //Final test for plausible decoding
    delta=60;                                //Fano threshold step

    fftwf_complex *fftin, *fftout;

    int32_t mettab[2][256];
    float bias=0.42;

    // setup metric table
    for(i=0; i<256; i++) {
        mettab[0][i]=round( 10*(metric_tables[2][i]-bias) );
        mettab[1][i]=round( 10*(metric_tables[2][255-i]-bias) );
    }

    FILE *fp_fftw_wisdom_file, *fhash;

    if ((fp_fftw_wisdom_file = fopen("wspr_wisdom.dat", "r"))) {  //Open FFTW wisdom
        fftwf_import_wisdom_from_file(fp_fftw_wisdom_file);
        fclose(fp_fftw_wisdom_file);
    }

    // Do windowed ffts over 2 symbols, stepped by half symbols
    int32_t nffts=4*floor(npoints/512)-1;
    fftin=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
    fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
    PLAN3 = fftwf_plan_dft_1d(512, fftin, fftout, FFTW_FORWARD, PATIENCE);

    float ps[512][nffts];
    float w[512];
    for(i=0; i<512; i++) {
        w[i]=sin(0.006147931*i);
    }

    if( options.usehashtable ) {
        char line[80], hcall[12];
        if( (fhash=fopen("hashtable.txt","r+")) ) {
            while (fgets(line, sizeof(line), fhash) != NULL) {
                sscanf(line,"%d %s",&nh,hcall);
                strcpy(hashtab+nh*13,hcall);
            }
        } else {
            fhash=fopen("hashtable.txt","w+");
        }
        fclose(fhash);
    }

    //*************** main loop starts here *****************
    for (ipass=0; ipass<options.npasses; ipass++) {

        if( ipass == 1 && uniques == 0 ) break;
        if( ipass == 1 ) {  //otherwise we bog down on the second pass
            options.quickmode = 1;
        }

        memset(ps,0.0, sizeof(float)*512*nffts);
        for (i=0; i<nffts; i++) {
            for(j=0; j<512; j++ ) {
                k=i*128+j;
                fftin[j][0]=idat[k] * w[j];
                fftin[j][1]=qdat[k] * w[j];
            }
            fftwf_execute(PLAN3);
            for (j=0; j<512; j++ ) {
                k=j+256;
                if( k>511 )
                    k=k-512;
                ps[j][i]=fftout[k][0]*fftout[k][0]+fftout[k][1]*fftout[k][1];
            }
        }

        sched_yield();

        // Compute average spectrum
        float psavg[512]= {0};
        for (i=0; i<nffts; i++) {
            for (j=0; j<512; j++) {
                psavg[j]=psavg[j]+ps[j][i];
            }
        }

        sched_yield();

        // Smooth with 7-point window and limit spectrum to +/-150 Hz
        int32_t window[7]= {1,1,1,1,1,1,1};
        float smspec[411];
        for (i=0; i<411; i++) {
            smspec[i]=0.0;
            for(j=-3; j<=3; j++) {
                k=256-205+i+j;
                smspec[i]=smspec[i]+window[j+3]*psavg[k];
            }
        }

        // Sort spectrum values, then pick off noise level as a percentile
        float tmpsort[411];
        for (j=0; j<411; j++) {
            tmpsort[j]=smspec[j];
        }
        qsort(tmpsort, 411, sizeof(float), floatcomp);

        // Noise level of spectrum is estimated as 123/411= 30'th percentile
        float noise_level = tmpsort[122];

        sched_yield();

        /* Renormalize spectrum so that (large) peaks represent an estimate of snr.
         * We know from experience that threshold snr is near -7dB in wspr bandwidth,
         * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
         * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15. */

        float min_snr, snr_scaling_factor;
        min_snr = pow(10.0,-7.0/10.0); //this is min snr in wspr bw
        if( wspr_type == 2 ) {
            snr_scaling_factor=26.3;
        } else {
            snr_scaling_factor=35.3;
        }
        for (j=0; j<411; j++) {
            smspec[j]=smspec[j]/noise_level - 1.0;
            if( smspec[j] < min_snr) smspec[j]=0.1;
            continue;
        }

        // Find all local maxima in smoothed spectrum.
        for (i=0; i<200; i++) {
            freq0[i]=0.0;
            snr0[i]=0.0;
            drift0[i]=0.0;
            shift0[i]=0;
            sync0[i]=0.0;
        }

        int32_t npk=0;
        for(j=1; j<410; j++) {
            if((smspec[j]>smspec[j-1]) && (smspec[j]>smspec[j+1]) && (npk<200)) {
                freq0[npk]=(j-205)*(DF/2.0);
                snr0[npk]=10*log10(smspec[j])-snr_scaling_factor;
                npk++;
            }
        }

        // Compute corrected fmin, fmax, accounting for dial frequency error
        fmin += dialfreq_error;    // dialfreq_error is in units of Hz
        fmax += dialfreq_error;

        // Don't waste time on signals outside of the range [fmin,fmax].
        i=0;
        for( j=0; j<npk; j++) {
            if( freq0[j] >= fmin && freq0[j] <= fmax ) {
                freq0[i]=freq0[j];
                snr0[i]=snr0[j];
                i++;
            }
        }
        npk=i;

        // bubble sort on snr, bringing freq along for the ride
        int32_t pass;
        float tmp;
        for (pass = 1; pass <= npk - 1; pass++) {
            for (k = 0; k < npk - pass ; k++) {
                if (snr0[k] < snr0[k+1]) {
                    tmp = snr0[k];
                    snr0[k] = snr0[k+1];
                    snr0[k+1] = tmp;
                    tmp = freq0[k];
                    freq0[k] = freq0[k+1];
                    freq0[k+1] = tmp;
                }
            }
        }

        sched_yield();

        /* Make coarse estimates of shift (DT), freq, and drift
         * Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
           to nominal start time, which is 2 seconds into the file
         * Calculates shift relative to the beginning of the file
         * Negative shifts mean that signal started before start of file
         * The program prints DT = shift-2 s
         * Shifts that cause sync vector to fall off of either end of the data
           vector are accommodated by "partial decoding", such that missing
           symbols produce a soft-decision symbol value of 128
         * The frequency drift model is linear, deviation of +/- drift/2 over the
           span of 162 symbols, with deviation equal to 0 at the center of the
           signal vector.
         */

        int32_t idrift,ifr,if0,ifd,k0;
        int32_t kindex;
        float smax,ss,pow,p0,p1,p2,p3;
        for(j=0; j<npk; j++) {                              //For each candidate...
            smax=-1e30;
            if0=freq0[j]/(DF/2.0)+NSPERSYM;
            for (ifr=if0-1; ifr<=if0+1; ifr++) {                      //Freq search
                for( k0=-10; k0<22; k0++) {                             //Time search
                    for (idrift=-maxdrift; idrift<=maxdrift; idrift++) {  //Drift search
                        ss=0.0;
                        pow=0.0;
                        for (k=0; k<NSYM; k++) {                             //Sum over symbols
                            ifd=ifr+((float)k-81.0)/81.0*( (float)idrift )/DF;
                            kindex=k0+2*k;
                            if( kindex < nffts ) {
                                p0=ps[ifd-3][kindex];
                                p1=ps[ifd-1][kindex];
                                p2=ps[ifd+1][kindex];
                                p3=ps[ifd+3][kindex];

                                p0=sqrtf(p0);
                                p1=sqrtf(p1);
                                p2=sqrtf(p2);
                                p3=sqrtf(p3);

                                ss=ss+(2*pr3[k]-1)*((p1+p3)-(p0+p2));
                                pow=pow+p0+p1+p2+p3;
                                sync1=ss/pow;
                            }
                        }
                        if( sync1 > smax ) {                  //Save coarse parameters
                            smax=sync1;
                            shift0[j]=128*(k0+1);
                            drift0[j]=idrift;
                            freq0[j]=(ifr-NSPERSYM)*(DF/2.0);
                            sync0[j]=sync1;
                        }
                    }
                }
            }
        }

        sched_yield();

        /*
         Refine the estimates of freq, shift using sync as a metric.
         Sync is calculated such that it is a float taking values in the range
         [0.0,1.0].

         Function sync_and_demodulate has three modes of operation
         mode is the last argument:

         0 = no frequency or drift search. find best time lag.
         1 = no time lag or drift search. find best frequency.
         2 = no frequency or time lag search. Calculate soft-decision
         symbols using passed frequency and shift.

         NB: best possibility for OpenMP may be here: several worker threads
         could each work on one candidate at a time.
         */

        for (j=0; j<npk; j++) {
            memset(symbols,0,sizeof(char)*NBITS*2);
            memset(callsign,0,sizeof(char)*13);
            memset(call_loc_pow,0,sizeof(char)*23);
            memset(call,0,sizeof(char)*13);
            memset(loc,0,sizeof(char)*7);
            memset(pwr,0,sizeof(char)*3);

            f1=freq0[j];
            drift1=drift0[j];
            shift1=shift0[j];
            sync1=sync0[j];

            // Fine search for best sync lag (mode 0)
            fstep=0.0;
            lagmin=shift1-144;
            lagmax=shift1+144;
            lagstep=8;

            if(options.quickmode)
                lagstep=16;

            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);

            // Fine search for frequency peak (mode 1)
            fstep=0.1;
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);

            if( sync1 > minsync1 ) {
                worth_a_try = 1;
            } else {
                worth_a_try = 0;
            }

            int32_t idt=0, ii=0, jiggered_shift;
            double y,sq,rms;
            not_decoded=1;

            while ( worth_a_try && not_decoded && idt<=(128/iifac)) {
                ii=(idt+1)/2;
                if( idt%2 == 1 ) ii=-ii;
                ii=iifac*ii;
                jiggered_shift=shift1+ii;

                // Use mode 2 to get soft-decision symbols
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep,
                                    &jiggered_shift, lagmin, lagmax, lagstep, &drift1, symfac,
                                    &sync1, 2);

                sq=0.0;
                for(i=0; i<NSYM; i++) {
                    y=(double)symbols[i] - 128.0;
                    sq += y*y;
                }
                rms=sqrt(sq/NSYM);

                if((sync1 > minsync2) && (rms > minrms)) {
                    deinterleave(symbols);
                    not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,NBITS,
                                       mettab,delta,maxcycles);
                }
                idt++;
                if( options.quickmode ) break;
            }

            if( worth_a_try && !not_decoded ) {
                for(i=0; i<11; i++) {
                    if( decdata[i]>127 ) {
                        message[i]=decdata[i]-256;
                    } else {
                        message[i]=decdata[i];
                    }
                }

                // Unpack the decoded message, update the hashtable, apply
                // sanity checks on grid and power, and return
                // call_loc_pow string and also callsign (for de-duping).
                noprint=unpk_(message,hashtab,call_loc_pow,call,loc,pwr,callsign);

                if( options.subtraction && (ipass == 0) && !noprint ) {
                    unsigned char channel_symbols[NSYM];

                    if( get_wspr_channel_symbols(call_loc_pow, hashtab, channel_symbols) ) {
                        subtract_signal2(idat, qdat, npoints, f1, shift1, drift1, channel_symbols);
                    } else {
                        break;
                    }

                }

                // Remove dupes (same callsign and freq within 3 Hz)
                int32_t dupe=0;
                for (i=0; i<uniques; i++) {
                    if(!strcmp(callsign,allcalls[i]) && (fabs(f1-allfreqs[i]) <3.0))
                        dupe=1;
                }

                if( (verbose || !dupe) && !noprint) {
                    strcpy(allcalls[uniques],callsign);
                    allfreqs[uniques]=f1;
                    uniques++;

                    // Add an extra space at the end of each line so that wspr-x doesn't
                    // truncate the power (TNX to DL8FCL!)
                    if( wspr_type == 15 ) {
                        freq_print=dialfreq+(1500+112.5+f1/8.0)/1e6;
                        dt_print=shift1*8*DT-2.0;
                    } else {
                        freq_print=dialfreq+(1500+f1)/1e6;
                        dt_print=shift1*DT-2.0;
                    }

                    decodes[uniques-1].sync=sync1;
                    decodes[uniques-1].snr=snr0[j];
                    decodes[uniques-1].dt=dt_print;
                    decodes[uniques-1].freq=freq_print;
                    decodes[uniques-1].drift=drift1;
                    decodes[uniques-1].cycles=cycles;
                    decodes[uniques-1].jitter=ii;
                    strcpy(decodes[uniques-1].message,call_loc_pow);
                    strcpy(decodes[uniques-1].call,call);
                    strcpy(decodes[uniques-1].loc,loc);
                    strcpy(decodes[uniques-1].pwr,pwr);
                }
            }
        }
		sched_yield();
    }

    // sort the result
    struct decoder_results temp;
    for (j = 1; j <= uniques - 1; j++) {
        for (k = 0; k < uniques - j ; k++) {
            if (decodes[k].snr < decodes[k+1].snr) {
                temp = decodes[k];
                decodes[k]=decodes[k+1];;
                decodes[k+1] = temp;
            }
        }
    }

    // Return number of spots to the calling fct
    *n_results = uniques;

    fftwf_free(fftin);
    fftwf_free(fftout);

    if ((fp_fftw_wisdom_file = fopen("wspr_wisdom.dat", "w"))) {
        fftwf_export_wisdom_to_file(fp_fftw_wisdom_file);
        fclose(fp_fftw_wisdom_file);
    }

    fftwf_destroy_plan(PLAN1);
    fftwf_destroy_plan(PLAN2);
    fftwf_destroy_plan(PLAN3);

    if( options.usehashtable ) {
        fhash=fopen("hashtable.txt","w");
        for (i=0; i<32768; i++) {
            if( strncmp(hashtab+i*13,"\0",1) != 0 ) {
                fprintf(fhash,"%5d %s\n",i,hashtab+i*13);
            }
        }
        fclose(fhash);
    }

    if(writenoise == 999) return -1;  //Silence compiler warning
    return 0;
}

/* Callback for each buffer received */
static void rtlsdr_callback(unsigned char *samples, uint32_t samples_count, void *ctx) {
    int8_t *sigIn = (int8_t*) samples;
    uint32_t sigLenght = samples_count;

    static uint32_t decimationIndex=0;

    /* CIC buffers */
    static int32_t  Ix1,Ix2,Qx1,Qx2;
    static int32_t  Iy1,It1y,It1z,Qy1,Qt1y,Qt1z;
    static int32_t  Iy2,It2y,It2z,Qy2,Qt2y,Qt2z;

    /* FIR compensation filter buffers */
    static float    firI[32], firQ[32];

    /* FIR compensation filter coefs
       Using : Octave/MATLAB code for generating compensation FIR coefficients
       URL : https://github.com/WestCoastDSP/CIC_Octave_Matlab
     */
    const static float zCoef[33] = {
        -0.0027772683, -0.0005058826,  0.0049745750, -0.0034059318,
        -0.0077557814,  0.0139375423,  0.0039896935, -0.0299394142,
         0.0162250643,  0.0405130860, -0.0580746013, -0.0272104968,
         0.1183705475, -0.0306029022, -0.2011241667,  0.1615898423,
         0.5000000000,
         0.1615898423, -0.2011241667, -0.0306029022,  0.1183705475,
        -0.0272104968, -0.0580746013,  0.0405130860,  0.0162250643,
        -0.0299394142,  0.0039896935,  0.0139375423, -0.0077557814,
        -0.0034059318,  0.0049745750, -0.0005058826, -0.0027772683
    };
    float Isum,Qsum;

    /* Convert unsigned to signed */
    for(uint32_t i=0; i<sigLenght; i++)
        sigIn[i] ^= 0x80;  // XOR with a binary mask to flip the first bit (sign)
        //sigIn[i] = (int8_t)((int32_t)samples[i] - 127);

    /* Economic mixer @ fs/4 (upper band)
       At fs/4, sin and cosin calculation are no longueur necessary.

               0   | pi/2 |  pi  | 3pi/2
             ----------------------------
       sin =   0   |  1   |  0   |  -1  |
       cos =   1   |  0   | -1   |   0  |

       out_I = in_I * cos(x) - in_Q * sin(x)
       out_Q = in_Q * cos(x) + in_I * sin(x)
       (Weaver technique, keep the upper band, IQ inverted on RTL devices)
    */
    int8_t tmp;
    for (uint32_t i=0; i<sigLenght; i+=8) {
        tmp = -sigIn[i+3];
        sigIn[i+3] = sigIn[i+2];
        sigIn[i+2] = tmp;

        sigIn[i+4] = -sigIn[i+4];
        sigIn[i+5] = -sigIn[i+5];

        tmp = -sigIn[i+6];
        sigIn[i+6] = sigIn[i+7];
        sigIn[i+7] = tmp;
    }

    /* CIC decimator (N=2)
       (could be not perfect in time for some sampling rate.
       Ex: AirSpy vs AirSpy Mini, but works fine in practice)
       Info: * Understanding CIC Compensation Filters
               https://www.altera.com/en_US/pdfs/literature/an/an455.pdf
             * Understanding cascaded integrator-comb filters
               http://www.embedded.com/design/configurable-systems/4006446/Understanding-cascaded-integrator-comb-filters
    */
    for(int32_t i=0; i<sigLenght/2; i++) {
        /* Integrator stages (N=2) */
        Ix1 += (int32_t)sigIn[i*2];
        Qx1 += (int32_t)sigIn[i*2+1];
        Ix2 += Ix1;
        Qx2 += Qx1;

        /* Decimation R=6400 */
        decimationIndex++;
        if (decimationIndex < DOWNSAMPLING) {
            continue;
        }

        // FIXME/TODO : some optimisition here
        /* 1st Comb */
        Iy1  = Ix2 - It1z;
        It1z = It1y;
        It1y = Ix2;
        Qy1  = Qx2 - Qt1z;
        Qt1z = Qt1y;
        Qt1y = Qx2;

        /* 2nd Comd */
        Iy2  = Iy1 - It2z;
        It2z = It2y;
        It2y = Iy1;
        Qy2  = Qy1 - Qt2z;
        Qt2z = Qt2y;
        Qt2y = Qy1;

        // FIXME/TODO : could be made with int32_t (8 bits, 20 bits)
        /* FIR compensation filter */
        Isum=0.0, Qsum=0.0;
        for (uint32_t j=0; j<32; j++) {
            Isum += firI[j]*zCoef[j];
            Qsum += firQ[j]*zCoef[j];
            if (j<31) {
                firI[j] = firI[j+1];
                firQ[j] = firQ[j+1];
            }
        }
        firI[31] = (float)Iy2;
        firQ[31] = (float)Qy2;
        Isum += firI[31]*zCoef[32];
        Qsum += firQ[31]*zCoef[32];

        /* Save the result in the buffer */
        if (rx_state.iqIndex < (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE)) {
            /* Lock the buffer during writing */     // Overkill ?!
            pthread_rwlock_wrlock(&dec.rw);
            rx_state.iSamples[rx_state.iqIndex] = Isum;
            rx_state.qSamples[rx_state.iqIndex] = Qsum;
            pthread_rwlock_unlock(&dec.rw);
            rx_state.iqIndex++;
        } else {
            if (rx_state.decode_flag == false) {
                /* Send a signal to the other thread to start the decoding */
                pthread_mutex_lock(&dec.ready_mutex);
                pthread_cond_signal(&dec.ready_cond);
                pthread_mutex_unlock(&dec.ready_mutex);
                rx_state.decode_flag = true;
                printf("RX done! [Buffer size: %d]\n", rx_state.iqIndex);
            }
        }
        decimationIndex = 0;
    }
}


/* Thread for RX blocking function */
static void *rtlsdr_rx(void *arg) {
    /* Read & blocking call */
    printf("rtlsdr Read & blocking call\n");
    rtlsdr_read_async(rtl_device, rtlsdr_callback, NULL, 0, DEFAULT_BUF_LENGTH);
    exit(0);
    return 0;
}


void postSpots(uint32_t n_results) {
    CURL *curl;
    CURLcode res;
    char url[256]; // FIXME, possible buffer overflow

    for (uint32_t i=0; i<n_results; i++) {
        sprintf(url,"http://wsprnet.org/post?function=wspr&rcall=%s&rgrid=%s&rqrg=%.6f&date=%s&time=%s&sig=%.0f&dt=%.1f&tqrg=%.6f&tcall=%s&tgrid=%s&dbm=%s&version=0.2r_wsprd&mode=2",
                dec_options.rcall, dec_options.rloc, dec_results[i].freq, dec_options.date, dec_options.uttime,
                dec_results[i].snr, dec_results[i].dt, dec_results[i].freq,
                dec_results[i].call, dec_results[i].loc, dec_results[i].pwr);

        printf("Spot : %3.2f %4.2f %10.6f %2d  %-s\n",
               dec_results[i].snr, dec_results[i].dt, dec_results[i].freq,
               (int)dec_results[i].drift, dec_results[i].message);

        curl = curl_easy_init();
        if(curl) {
            curl_easy_setopt(curl, CURLOPT_URL, url);
            curl_easy_setopt(curl, CURLOPT_NOBODY, 1);
            res = curl_easy_perform(curl);

            if(res != CURLE_OK)
                fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));

            curl_easy_cleanup(curl);
        }
    }
    if (n_results == 0)
        printf("No spot\n");
}


static void *wsprDecoder(void *arg) {
    /* WSPR decoder use buffers of 45000 samples (hardcoded)
       (120 sec max @ 375sps = 45000 samples)
    */
    static float iSamples[45000]={0};
    static float qSamples[45000]={0};
    static uint32_t samples_len;
    int32_t n_results=0;

    while (!rx_state.exit_flag) {
        pthread_mutex_lock(&dec.ready_mutex);
        pthread_cond_wait(&dec.ready_cond, &dec.ready_mutex);
        pthread_mutex_unlock(&dec.ready_mutex);

        if(rx_state.exit_flag)  // Abord case, final sig
            break;

        /* Lock the buffer access and make a local copy */
        pthread_rwlock_wrlock(&dec.rw);
        memcpy(iSamples, rx_state.iSamples, rx_state.iqIndex * sizeof(float));
        memcpy(qSamples, rx_state.qSamples, rx_state.iqIndex * sizeof(float));
        samples_len = rx_state.iqIndex;  // Overkill ?
        pthread_rwlock_unlock(&dec.rw);

        /* Date and time will be updated/overload during the search & decoding process
           Make a simple copy
        */
        memcpy(dec_options.date, rx_options.date, sizeof(rx_options.date));
        memcpy(dec_options.uttime, rx_options.uttime, sizeof(rx_options.uttime));

        /* DEBUG -- Save samples
        printf("Writing file\n");
        FILE* fd = NULL;
        fd = fopen("samples.bin", "wb");
        int r=fwrite(rx_state.iSamples, sizeof(float), samples_len, fd);
        printf("%d samples written file\n", r);
        fclose(fd);
        */

        /* Search & decode the signal */
        printf("Search & decode the signal\n");
        wspr_decode(iSamples, qSamples, samples_len, dec_options, dec_results, &n_results);
        printf("Going to post spot\n");
        postSpots(n_results);

    }
    pthread_exit(NULL);
}


double atofs(char *s) {
    /* standard suffixes */
    char last;
    uint32_t len;
    double suff = 1.0;
    len = strlen(s);
    last = s[len-1];
    s[len-1] = '\0';
    switch (last) {
    case 'g':
    case 'G':
        suff *= 1e3;
    case 'm':
    case 'M':
        suff *= 1e3;
    case 'k':
    case 'K':
        suff *= 1e3;
        suff *= atof(s);
        s[len-1] = last;
        return suff;
    }
    s[len-1] = last;
    return atof(s);
}


int32_t parse_u64(char* s, uint64_t* const value) {
    uint_fast8_t base = 10;
    char* s_end;
    uint64_t u64_value;

    if( strlen(s) > 2 ) {
        if( s[0] == '0' ) {
            if( (s[1] == 'x') || (s[1] == 'X') ) {
                base = 16;
                s += 2;
            } else if( (s[1] == 'b') || (s[1] == 'B') ) {
                base = 2;
                s += 2;
            }
        }
    }

    s_end = s;
    u64_value = strtoull(s, &s_end, base);
    if( (s != s_end) && (*s_end == 0) ) {
        *value = u64_value;
        return 1;
    } else {
        return 0;
    }
}


/* Reset flow control variable & decimation variables */
void initSampleStorage() {
    rx_state.decode_flag = false;
    rx_state.iqIndex=0;
}


/* Default options for the decoder */
void initDecoder_options() {
    dec_options.usehashtable = 1;
    dec_options.npasses = 2;
    dec_options.subtraction = 1;
    dec_options.quickmode = 0;
}


/* Default options for the receiver */
void initrx_options() {
    rx_options.gain = 29;
    rx_options.autogain = 0;
    rx_options.ppm = 0;
    rx_options.shift = 0;
    rx_options.directsampling = 0;
}


void sigint_callback_handler(int signum) {
    fprintf(stdout, "Caught signal %d\n", signum);
    rx_state.exit_flag = true;
}


void usage(void) {
    fprintf(stderr,
            "rtlsdr_wsprd, a simple WSPR daemon for RTL receivers\n\n"
            "Use:\trtlsdr_wsprd -f frequency -c callsign -l locator [options]\n"
            "\t-f dial frequency [(,k,M) Hz], check http://wsprnet.org/ for freq.\n"
            "\t-c your callsign (12 chars max)\n"
            "\t-l your locator grid (6 chars max)\n"
            "Receiver extra options:\n"
            "\t-g gain [0-49] (default: 29)\n"
            "\t-a auto gain (default: off)\n"
            "\t-o frequency offset (default: 0)\n"
            "\t-p crystal correction factor (ppm) (default: 0)\n"
            "\t-u upconverter (default: 0, example: 125M)\n"
            "\t-d direct dampling [0,1,2] (default: 0, 1 for I input, 2 for Q input)\n"
            "Decoder extra options:\n"
            "\t-H do not use (or update) the hash table\n"
            "\t-Q quick mode, doesn't dig deep for weak signals\n"
            "\t-S single pass mode, no subtraction (same as original wsprd)\n"
            "Example:\n"
            "\trtlsdr_wsprd -f 144.489M -c A1XYZ -l AB12cd -g 29 -o -4200\n");
    exit(1);
}


int main(int argc, char **argv) {
    uint32_t opt;

    int32_t  rtl_result;
    uint32_t rtl_index = 0; // By default, use the first RTLSDR
    int32_t  rtl_count;
    char     rtl_vendor[256], rtl_product[256], rtl_serial[256];

    initrx_options();
    initDecoder_options();

    /* RX buffer allocation */
    rx_state.iSamples=malloc(sizeof(float)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);
    rx_state.qSamples=malloc(sizeof(float)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);

    /* Stop condition setup */
    rx_state.exit_flag   = false;
    rx_state.decode_flag = false;

    if (argc <= 1)
        usage();

    while ((opt = getopt(argc, argv, "f:c:l:g:a:o:p:u:d:H:Q:S")) != -1) {
        switch (opt) {
        case 'f': // Frequency
            rx_options.dialfreq = (uint32_t)atofs(optarg);
            break;
        case 'c': // Callsign
            sprintf(dec_options.rcall, "%.12s", optarg);
            break;
        case 'l': // Locator / Grid
            sprintf(dec_options.rloc, "%.6s", optarg);
            break;
        case 'g': // Small signal amplifier gain
            rx_options.gain = atoi(optarg);
            if (rx_options.gain < 0) rx_options.gain = 0;
            if (rx_options.gain > 49 ) rx_options.gain = 49;
            rx_options.gain *= 10;
            break;
        case 'a': // Auto gain
            rx_options.autogain = atoi(optarg);
            if (rx_options.autogain < 0) rx_options.autogain = 0;
            if (rx_options.autogain > 1 ) rx_options.autogain = 1;
            break;
        case 'o': // Fine frequency correction
            rx_options.shift = atoi(optarg);
            break;
        case 'p':
            rx_options.ppm = atoi(optarg);
            break;
        case 'u': // Upconverter frequency
            rx_options.upconverter = (uint32_t)atofs(optarg);
            break;
        case 'd': // Direct Sampling
            rx_options.directsampling = (uint32_t)atofs(optarg);
            break;
        case 'H': // Decoder option, use a hastable
            dec_options.usehashtable = 0;
            break;
        case 'Q': // Decoder option, faster
            dec_options.quickmode = 1;
            break;
        case 'S': // Decoder option, single pass mode (same as original wsprd)
            dec_options.subtraction = 0;
            dec_options.npasses = 1;
            break;
        default:
            usage();
            break;
        }
    }

    if (rx_options.dialfreq == 0) {
        fprintf(stderr, "Please specify a dial frequency.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    if (dec_options.rcall[0] == 0) {
        fprintf(stderr, "Please specify your callsign.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    if (dec_options.rloc[0] == 0) {
        fprintf(stderr, "Please specify your locator.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    /* Calcule shift offset */
    rx_options.realfreq = rx_options.dialfreq + rx_options.shift + rx_options.upconverter;

    /* Store the frequency used for the decoder */
    dec_options.freq = rx_options.dialfreq;

    /* If something goes wrong... */
    signal(SIGINT, &sigint_callback_handler);
    signal(SIGTERM, &sigint_callback_handler);
    signal(SIGILL, &sigint_callback_handler);
    signal(SIGFPE, &sigint_callback_handler);
    signal(SIGSEGV, &sigint_callback_handler);
    signal(SIGABRT, &sigint_callback_handler);

    /* Init & parameter the device */
    rtl_count = rtlsdr_get_device_count();
    if (!rtl_count) {
        fprintf(stderr, "No supported devices found\n");
        return EXIT_FAILURE;
    }


    fprintf(stderr, "Found %d device(s):\n", rtl_count);
    for (uint32_t i=0; i<rtl_count; i++) {
        rtlsdr_get_device_usb_strings(i, rtl_vendor, rtl_product, rtl_serial);
        fprintf(stderr, "  %d:  %s, %s, SN: %s\n", i, rtl_vendor, rtl_product, rtl_serial);
    }
    fprintf(stderr, "\nUsing device %d: %s\n", rtl_index, rtlsdr_get_device_name(rtl_index));


    rtl_result = rtlsdr_open(&rtl_device, rtl_index);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to open rtlsdr device #%d.\n", rtl_index);
        return EXIT_FAILURE;
    }

    if (rx_options.directsampling) {
        rtl_result = rtlsdr_set_direct_sampling(rtl_device, rx_options.directsampling);
        if (rtl_result < 0) {
            fprintf(stderr, "ERROR: Failed to set direct sampling\n");
            rtlsdr_close(rtl_device);
            return EXIT_FAILURE;
        }
    }

    rtl_result = rtlsdr_set_sample_rate(rtl_device, SAMPLING_RATE);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to set sample rate\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }


    rtl_result = rtlsdr_set_tuner_gain_mode(rtl_device, 1);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to enable manual gain\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }


    if (rx_options.autogain) {
        rtl_result = rtlsdr_set_tuner_gain_mode(rtl_device, 0);
        if (rtl_result != 0) {
            fprintf(stderr, "ERROR: Failed to set tuner gain\n");
            rtlsdr_close(rtl_device);
            return EXIT_FAILURE;
        }
    } else {
        rtl_result = rtlsdr_set_tuner_gain(rtl_device, rx_options.gain);
        if (rtl_result != 0) {
            fprintf(stderr, "ERROR: Failed to set tuner gain\n");
            rtlsdr_close(rtl_device);
            return EXIT_FAILURE;
        }
    }


    if (rx_options.ppm != 0) {
        rtl_result = rtlsdr_set_freq_correction(rtl_device, rx_options.ppm);
        if (rtl_result < 0) {
            fprintf(stderr, "ERROR: Failed to set ppm error\n");
            rtlsdr_close(rtl_device);
            return EXIT_FAILURE;
        }
    }


    rtl_result = rtlsdr_set_center_freq(rtl_device, rx_options.realfreq + FS4_RATE + 1500);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to set frequency\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }


    rtl_result = rtlsdr_reset_buffer(rtl_device);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to reset buffers.\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }

    /* Print used parameter */
    time_t rawtime;
    time ( &rawtime );
    struct tm *gtm = gmtime(&rawtime);
    printf("\nStarting rtlsdr-wsprd (%04d-%02d-%02d, %02d:%02dz) -- Version 0.2\n",
           gtm->tm_year + 1900, gtm->tm_mon + 1, gtm->tm_mday, gtm->tm_hour, gtm->tm_min);
    printf("  Callsign     : %s\n", dec_options.rcall);
    printf("  Locator      : %s\n", dec_options.rloc);
    printf("  Dial freq.   : %d Hz\n", rx_options.dialfreq);
    printf("  Real freq.   : %d Hz\n", rx_options.realfreq);
    printf("  PPM factor   : %d\n", rx_options.ppm);
    if(rx_options.autogain)
        printf("  Auto gain    : enable\n");
    else
        printf("  Gain         : %d dB\n", rx_options.gain/10);


    /* Time alignment stuff */
    struct timeval lTime;
    gettimeofday(&lTime, NULL);
    uint32_t sec   = lTime.tv_sec % 120;
    uint32_t usec  = sec * 1000000 + lTime.tv_usec;
    uint32_t uwait = 120000000 - usec;
    printf("Wait for time sync (start in %d sec)\n\n", uwait/1000000);

    /* Prepare a low priority param for the decoder thread */
    struct sched_param param;
    pthread_attr_init(&dec.tattr);
    pthread_attr_setschedpolicy(&dec.tattr, SCHED_RR);
    pthread_attr_getschedparam(&dec.tattr, &param);
    param.sched_priority = 90;  // = sched_get_priority_min();
    pthread_attr_setschedparam(&dec.tattr, &param);
    //int res=0;
    //printf("get: %d\n", res)

    /* Create a thread and stuff for separate decoding
       Info : https://computing.llnl.gov/tutorials/pthreads/
    */
    pthread_rwlock_init(&dec.rw, NULL);
    pthread_cond_init(&dec.ready_cond, NULL);
    pthread_mutex_init(&dec.ready_mutex, NULL);
    pthread_create(&dongle.thread, NULL, rtlsdr_rx, NULL);
    pthread_create(&dec.thread, &dec.tattr, wsprDecoder, NULL);

    /* Main loop : Wait, read, decode */
    while (!rx_state.exit_flag) {
        /* Wait for time Sync on 2 mins */
        gettimeofday(&lTime, NULL);
        sec   = lTime.tv_sec % 120;
        usec  = sec * 1000000 + lTime.tv_usec;
        uwait = 120000000 - usec + 10000;  // Adding 10ms, to be sure to reach this next minute
        usleep(uwait);
        printf("SYNC! RX started\n");

        /* Use the Store the date at the begin of the frame */
        time ( &rawtime );
        gtm = gmtime(&rawtime);
        sprintf(rx_options.date,"%02d%02d%02d", gtm->tm_year - 100, gtm->tm_mon + 1, gtm->tm_mday);
        sprintf(rx_options.uttime,"%02d%02d", gtm->tm_hour, gtm->tm_min);

        /* Start to store the samples */
        printf("Start to store the samples\n");
        initSampleStorage();

        while( (rx_state.exit_flag == false) &&
               (rx_state.iqIndex < (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE) ) ) {
            usleep(250000);
        }
    }

    /* Stop the RX and free the blocking function */
    rtlsdr_cancel_async(rtl_device);

    /* Close the RTL device */
    rtlsdr_close(rtl_device);

    printf("Bye!\n");

    /* Wait the thread join (send a signal before to terminate the job) */
    pthread_mutex_lock(&dec.ready_mutex);
    pthread_cond_signal(&dec.ready_cond);
    pthread_mutex_unlock(&dec.ready_mutex);
    pthread_join(dec.thread, NULL);
    pthread_join(dongle.thread, NULL);

    /* Destroy the lock/cond/thread */
    pthread_rwlock_destroy(&dec.rw);
    pthread_cond_destroy(&dec.ready_cond);
    pthread_mutex_destroy(&dec.ready_mutex);
    pthread_exit(NULL);

    return EXIT_SUCCESS;
}
