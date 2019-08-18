
#ifndef ADBUFFER_LOADED
#define ADBUFFER_LOADED 1

extern void printtraffic() ;

extern void pushinteger4_(int x)  ;
extern void lookinteger4_(int *x) ;
extern void popinteger4_(int *x) ;

extern void pushcontrol1b_(int cc) ;
extern void lookcontrol1b_(int *cc) ;
extern void popcontrol1b_(int *cc) ;
extern void pushcontrol2b_(int cc) ;
extern void lookcontrol2b_(int *cc) ;
extern void popcontrol2b_(int *cc) ;
extern void pushcontrol3b_(int cc) ;
extern void lookcontrol3b_(int *cc) ;
extern void popcontrol3b_(int *cc) ;
extern void pushcontrol4b_(int cc) ;
extern void lookcontrol4b_(int *cc) ;
extern void popcontrol4b_(int *cc) ;
extern void pushcontrol5b_(int cc) ;
extern void lookcontrol5b_(int *cc) ;
extern void popcontrol5b_(int *cc) ;

extern void pushreal4_(float x) ;
extern void lookreal4_(float *x) ;
extern void popreal4_(float *x) ;

extern void pushreal8_(float x) ;
extern void lookreal8_(float *x) ;
extern void popreal8_(float *x) ;

extern void printbuffertop() ;
extern void showallstacks_() ;

#endif
