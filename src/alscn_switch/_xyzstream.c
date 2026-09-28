#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { FILE *fp; long frame; } Stream;
static void stream_free(PyObject *cap) {
    Stream *s = (Stream *)PyCapsule_GetPointer(cap, "alscn_switch.stream");
    if (s) { if (s->fp) fclose(s->fp); free(s); }
    else PyErr_Clear();
}
static char *trim(char *s) { while (isspace((unsigned char)*s)) s++; char *e=s+strlen(s); while(e>s&&isspace((unsigned char)e[-1])) *--e=0; return s; }
static int prop_index(char *header, const char *name, int *width) {
    char *p=strstr(header,"Properties="); if(!p) return -1; p+=11; char quote=(*p=='\"'||*p=='\'')?*p++:0;
    char *end=p; while(*end && (!quote || *end!=quote)) end++; size_t n=(size_t)(end-p); char *spec=(char*)malloc(n+1); if(!spec)return -1; memcpy(spec,p,n);spec[n]=0;
    int column=0,result=-1; char *save=NULL; char *tok=strtok_r(spec,":",&save); while(tok){char *nm=tok;char *ty=strtok_r(NULL,":",&save);char *cnt=strtok_r(NULL,":",&save);if(!ty||!cnt)break;int w=atoi(cnt);if(!strcmp(nm,name)){result=column;*width=w;break;}column+=w;tok=strtok_r(NULL,":",&save);}
    free(spec);return result;
}
static PyObject *open_stream(PyObject *self, PyObject *arg) {
    const char *path; if(!PyArg_ParseTuple(arg,"s",&path))return NULL; FILE *fp=fopen(path,"r");if(!fp)return PyErr_SetFromErrnoWithFilename(PyExc_OSError,path);
    Stream *s=(Stream*)calloc(1,sizeof(Stream));if(!s){fclose(fp);return PyErr_NoMemory();}s->fp=fp;
    return PyCapsule_New(s,"alscn_switch.stream",stream_free);
}
static PyObject *next_frame(PyObject *self, PyObject *arg) {
    PyObject *cap; if(!PyArg_ParseTuple(arg,"O",&cap))return NULL;Stream *s=(Stream*)PyCapsule_GetPointer(cap,"alscn_switch.stream");if(!s)return NULL;
    char *line=NULL;size_t capline=0;ssize_t got;do{got=getline(&line,&capline,s->fp);if(got<0){free(line);Py_RETURN_NONE;}}while(got==1&&line[0]=='\n');
    char *ep=NULL;long nl=strtol(line,&ep,10);if(ep==line||nl<1||nl>10000000){free(line);PyErr_SetString(PyExc_ValueError,"invalid extxyz atom count");return NULL;}int n=(int)nl;
    got=getline(&line,&capline,s->fp);if(got<0){free(line);PyErr_SetString(PyExc_ValueError,"truncated extxyz header");return NULL;}
    char *header=strdup(line);if(!header){free(line);return PyErr_NoMemory();}
    int w=0,pos=prop_index(header,"pos",&w);if(pos<0||w!=3){free(header);free(line);PyErr_SetString(PyExc_ValueError,"extxyz must declare Properties with pos:R:3");return NULL;}
    int uw=0,unwrapped=prop_index(header,"unwrapped_position",&uw);if(uw!=3)unwrapped=-1;
    int cw=0,charge=prop_index(header,"charge",&cw);if(cw!=1)charge=-1;
    int fw=0,forces=prop_index(header,"forces",&fw);if(fw!=3)forces=-1;
    int bw=0,bec=prop_index(header,"bec",&bw);if(bw!=9)bec=-1;
    int speciesw=0,species=prop_index(header,"species",&speciesw);if(species<0)species=prop_index(header,"element",&speciesw);if(species<0||speciesw!=1){free(header);free(line);PyErr_SetString(PyExc_ValueError,"extxyz must declare species:S:1");return NULL;}
    double lattice[9]={0},time=NAN;char *lp=strstr(header,"Lattice=");if(lp){lp+=8;if(*lp=='\"'||*lp=='\'')lp++;for(int i=0;i<9;i++){lattice[i]=strtod(lp,&lp);}}
    char *tp=strstr(header,"Time=");if(tp)time=strtod(tp+5,NULL);
    char *sp=(char*)malloc((size_t)n);double *xyz=(double*)malloc((size_t)n*3*sizeof(double));double *uwxyz=(double*)malloc((size_t)n*3*sizeof(double));double *q=(double*)malloc((size_t)n*sizeof(double));double *f=(double*)malloc((size_t)n*3*sizeof(double));double *be=(double*)malloc((size_t)n*9*sizeof(double));
    if(!sp||!xyz||!uwxyz||!q||!f||!be){free(header);free(line);free(sp);free(xyz);free(uwxyz);free(q);free(f);free(be);return PyErr_NoMemory();}
    for(int i=0;i<n;i++){got=getline(&line,&capline,s->fp);if(got<0){PyErr_SetString(PyExc_ValueError,"truncated extxyz atom data");goto fail;}
        char *save=NULL;int col=0,coord=0;char *tok=strtok_r(line," \t\r\n",&save);char elem[16]="";while(tok){if(col==species)snprintf(elem,sizeof elem,"%s",tok);if(col>=pos&&col<pos+3)xyz[i*3+coord++]=strtod(tok,NULL);if(col>=unwrapped&&unwrapped>=0&&col<unwrapped+3)uwxyz[i*3+(col-unwrapped)]=strtod(tok,NULL);if(col==charge&&charge>=0)q[i]=strtod(tok,NULL);if(col>=forces&&forces>=0&&col<forces+3)f[i*3+(col-forces)]=strtod(tok,NULL);if(col>=bec&&bec>=0&&col<bec+9)be[i*9+(col-bec)]=strtod(tok,NULL);col++;tok=strtok_r(NULL," \t\r\n",&save);}
        if(!elem[0]||coord!=3){PyErr_SetString(PyExc_ValueError,"malformed atom row or missing position columns");goto fail;}
        sp[i]=!strcmp(elem,"Al")?1:!strcmp(elem,"Sc")?2:!strcmp(elem,"N")?3:0;
        if(unwrapped<0)memcpy(uwxyz+i*3,xyz+i*3,3*sizeof(double));if(charge<0)q[i]=NAN;if(forces<0)for(int k=0;k<3;k++)f[i*3+k]=NAN;if(bec<0)for(int k=0;k<9;k++)be[i*9+k]=NAN;
    }
    PyObject *d=PyDict_New();if(!d)goto fail;
    PyObject *v=PyBytes_FromStringAndSize(sp,n);PyDict_SetItemString(d,"species",v);Py_DECREF(v);
    v=PyBytes_FromStringAndSize((char*)xyz,(Py_ssize_t)n*3*sizeof(double));PyDict_SetItemString(d,"positions",v);Py_DECREF(v);
    v=PyBytes_FromStringAndSize((char*)uwxyz,(Py_ssize_t)n*3*sizeof(double));PyDict_SetItemString(d,"unwrapped_positions",v);Py_DECREF(v);
    v=PyBytes_FromStringAndSize((char*)q,(Py_ssize_t)n*sizeof(double));PyDict_SetItemString(d,"charges",v);Py_DECREF(v);
    v=PyBytes_FromStringAndSize((char*)f,(Py_ssize_t)n*3*sizeof(double));PyDict_SetItemString(d,"forces",v);Py_DECREF(v);
    v=PyBytes_FromStringAndSize((char*)be,(Py_ssize_t)n*9*sizeof(double));PyDict_SetItemString(d,"bec",v);Py_DECREF(v);
    v=Py_BuildValue("(ddddddddd)",lattice[0],lattice[1],lattice[2],lattice[3],lattice[4],lattice[5],lattice[6],lattice[7],lattice[8]);PyDict_SetItemString(d,"lattice",v);Py_DECREF(v);
    v=Py_BuildValue("d",time);PyDict_SetItemString(d,"time_fs",v);Py_DECREF(v);
    v=PyLong_FromLong(n);PyDict_SetItemString(d,"natoms",v);Py_DECREF(v);v=PyLong_FromLong(++s->frame);PyDict_SetItemString(d,"frame",v);Py_DECREF(v);
    free(header);free(line);free(sp);free(xyz);free(uwxyz);free(q);free(f);free(be);return d;
fail: free(header);free(line);free(sp);free(xyz);free(uwxyz);free(q);free(f);free(be);return NULL;
}
static PyObject *close_stream(PyObject *self, PyObject *arg) {
    PyObject *cap; if(!PyArg_ParseTuple(arg,"O",&cap))return NULL;
    Stream *s=(Stream*)PyCapsule_GetPointer(cap,"alscn_switch.stream");if(!s)return NULL;
    if(s->fp){fclose(s->fp);s->fp=NULL;}Py_RETURN_NONE;
}
static PyMethodDef methods[]={{"open_stream",open_stream,METH_VARARGS,"Open extxyz stream."},{"next_frame",next_frame,METH_VARARGS,"Read one extxyz frame."},{"close_stream",close_stream,METH_VARARGS,"Close extxyz stream."},{NULL,NULL,0,NULL}};
static struct PyModuleDef module={PyModuleDef_HEAD_INIT,"_xyzstream",NULL,-1,methods};
PyMODINIT_FUNC PyInit__xyzstream(void){return PyModule_Create(&module);}
