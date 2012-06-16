#ifndef PARI_ERRORS_H
#define PARI_ERRORS_H

/* global flag set whenever setjmp is called */
int setjmp_active;

/* global variable which holds the current pari error number */
int pari_error_number;

/* counter to check whether our sig-ons balance our sig-offs. */
int sig_on_sig_off;
 
/* declarations of PARI's error callbacks */
void (*cb_pari_ask_confirm)(const char *);
int  (*cb_pari_handle_exception)(long);
int  (*cb_pari_whatnow)(PariOUT *out, const char *, int);
void (*cb_pari_sigint)(void);
void (*cb_pari_err_recover)(long);

void set_error_handler( int (*handler)(long) ) {
  cb_pari_handle_exception = handler;
}

void set_error_recoverer( void (*recoverer)(long) ) {
  cb_pari_err_recover = recoverer;
}

#define SIG_ON_MACRO() {			\
    sig_on_sig_off += 1;			\
    setjmp_active = 1;				\
    if ( setjmp(jmp_env) ) {			\
      return NULL;				\
    }						\
  }						\

#endif
