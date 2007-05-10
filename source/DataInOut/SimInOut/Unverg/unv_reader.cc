#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>


#include "unv_dat.hh"

int ReadNodesUnivFile(int **nodeReverseLabels, Widget w, unv *unvf) {

  int   i,j,status,maxNodeLabel=0;
  char  buffer[255];
  int   topRange = MAX_NUM_NODES;
  dataset_781 set781;

  NODE_LABELS  = (int*) XtRealloc((char*) NODE_LABELS, 
                                  MAX_NUM_NODES*sizeof(int*));
  MAX_NODE_LABEL = 0;
    
  set781.connect(unvf);
  sprintf(buffer, "Reading nodal coordinates:");
  wprint(work,buffer);
  do {
    if (!set781.read_data(&status)) {
      ErrorDialog(XtParent(w),
                  "An error occured reading nodal coordinate data\n"
                  "Check universal file!");
      TimeoutCursors(FALSE, FALSE);
      return 0;
    }
    N_NODES++;
    if (N_NODES>=topRange) {
      topRange += MAX_NUM_NODES;
      NODES = (point*) XtRealloc(
        (char*) NODES, topRange*sizeof(point));
      NODE_LABELS = (int*) XtRealloc( (char*) NODE_LABELS,
                                      topRange*sizeof(int));
    }
    NODES[N_NODES].x1=(double) set781.data.x1;
    NODES[N_NODES].x2=(double) set781.data.x2;
    NODES[N_NODES].x3=(double) set781.data.x3;
    NODE_LABELS[N_NODES] = set781.data.node_label;
    if (set781.data.node_label>= maxNodeLabel)
      maxNodeLabel = set781.data.node_label;
    if ((N_NODES % 500)==0) {
      CheckForInterrupt();
    }
  } while (status!=dataset::SETEmpty);

  sprintf(buffer, " done.\n");
  wprint(work,buffer);
  sprintf(buffer, " Read total of %d nodes.\n",N_NODES);
  wprint(work,buffer);

  CheckForInterrupt();

  /* establish internal->external node number table */

  *nodeReverseLabels = (int*) XtMalloc(
    (maxNodeLabel+1)*sizeof(int));
  for (i=1; i<=maxNodeLabel; i++)
    (*nodeReverseLabels)[i] = 0;

  for (i=1; i<=N_NODES;i++) {
    j = NODE_LABELS[i];
    (*nodeReverseLabels)[j] = i;
  }

  //XtFree((char*) nodeLabels);

  NODES_X_MIN = NODES_Y_MIN =  1.e10;
  NODES_X_MAX = NODES_Y_MAX = -1.e10;
  for (i=1; i<=N_NODES; ++i) {
    if (NODES_X_MIN>NODES[i].x1) NODES_X_MIN=NODES[i].x1;
    if (NODES_X_MAX<NODES[i].x1) NODES_X_MAX=NODES[i].x1;
    if (NODES_Y_MIN>NODES[i].x2) NODES_Y_MIN=NODES[i].x2;
    if (NODES_Y_MAX<NODES[i].x2) NODES_Y_MAX=NODES[i].x2;
  }

  if (!set781.close_set()) {
    ErrorDialog(XtParent(w),
                "An error occured terminating nodal data input!\n"
                "Check universal file!");
    TimeoutCursors(FALSE, FALSE);
    return 0;
  }
  return maxNodeLabel;
}

int ReadElementsUnivFile(int *nodeReverseLabels, Widget w, unv *unvf, 
                         int maxNodeLabel, int **elemReverseLabels) {

  int nodesNotFound=0,status,i,j,maxElemLabel=0;
  int topRange = MAX_NUM_ELEMENTS;
  dataset_780 set780;
  char        buffer[256];

  set780.connect(unvf);
  sprintf(buffer, "Reading element data:");
  wprint(work,buffer);
  do {
    if (!set780.read_data(&status)) {
      if (status!=dataset::SETEmptySet) {
        ErrorDialog(XtParent(w),
                    "An errror occured reading element data!\n"
                    " Check universal file!");
        TimeoutCursors(FALSE, FALSE);
        return 0;
      }
      else {
        ErrorDialog(XtParent(w),
                    "Empty element topology data set found!\n"
                    " Continue with processing of universal file!");
        status = dataset::SETEmpty;
      }
    }
    else {
      N_ELEMENTS++;
      if (N_ELEMENTS>=topRange) {
        topRange += MAX_NUM_ELEMENTS;
        ELEMENTS = (element_data*) XtRealloc(
          (char*) ELEMENTS, topRange*sizeof(element_data));
      }
      ELEMENTS[N_ELEMENTS].label   = set780.data.element_label;
      ELEMENTS[N_ELEMENTS].color   = set780.data.color;

      int triangle = 0, nn_nodes = set780.data.n_nodes;
      /*
        if (set780.data.n_nodes==3) {
        triangle = 1;
        nn_nodes = 4;
        }
        if (set780.data.n_nodes==6) {
        triangle = 2;
        nn_nodes = 8;
        }
      */

      ELEMENTS[N_ELEMENTS].n_nodes = nn_nodes;
      ELEMENTS[N_ELEMENTS].points  = (int*)
                                     XtMalloc(nn_nodes*sizeof(int));
      ELEMENTS[N_ELEMENTS].fe_type = set780.data.fe_desc_id;
      //ELEMENTS[N_ELEMENTS].edges   = (int*) XtMalloc(nn_nodes*sizeof(int));
      for (i=0; i<set780.data.n_nodes; ++i) {
        if (set780.data.node_labels[i]<=maxNodeLabel) {
          j = nodeReverseLabels[set780.data.node_labels[i]];
          if (j==0) 
            nodesNotFound++;
          else
            ELEMENTS[N_ELEMENTS].points[i] = j;
        }
        else
          nodesNotFound++;
      }

      /*
        if (triangle==1) {
        ELEMENTS[N_ELEMENTS].points[3] = 
        ELEMENTS[N_ELEMENTS].points[0] ;
        }
        if (triangle==2) {
        ELEMENTS[N_ELEMENTS].points[7] = 
        ELEMENTS[N_ELEMENTS].points[0] ;
        ELEMENTS[N_ELEMENTS].points[6] = 
        ELEMENTS[N_ELEMENTS].points[0] ;
        }
      */
      if (set780.data.element_label>= maxElemLabel)
        maxElemLabel = set780.data.element_label;
      if ((N_ELEMENTS % 500)==0) {
        CheckForInterrupt();
      }
    }
  } while (status!=dataset::SETEmpty);

  sprintf(buffer, " done.\n");
  wprint(work,buffer);
  sprintf(buffer, " Read total of %d elements.\n",N_ELEMENTS);
  wprint(work,buffer);

  if (nodesNotFound>0) {
    sprintf(buffer,"%s %d %s",
            " ERROR: element definitions referenced\n",
            nodesNotFound, " undefined nodes\n");
    wprint(work,buffer);
    TimeoutCursors(FALSE, FALSE);
    return 0;
  }

  /* establish internal->external element number table */

  *elemReverseLabels = (int*) XtMalloc(
    (maxElemLabel+1)*sizeof(int));
  for (i=1; i<=maxElemLabel; i++)
    (*elemReverseLabels)[i] = 0;

  for (i=1; i<=N_ELEMENTS;i++) {
    j = ELEMENTS[i].label;
    (*elemReverseLabels)[j] = i;
  }

  CheckForInterrupt();
  if (!set780.close_set()) {
    ErrorDialog(XtParent(w),
                "An error occured terminating element data input!\n"
                "Check universal file!");
    return 0;
  }
  return 1;
}

int ReadNResultsUnivFile(int &i, Widget w, unv *unvf, int &blocks,
                         int *nodeReverseLabels) {

  dataset_55 set55;
  char       buffer[256];
  int        n,type_size,j,k=0,stopReading=0,status;

  set55.connect(unvf);
  ++i;
  sprintf(buffer, "Reading nodal results no. %d:",i+1);
  wprint(work,buffer);
  if (!set55.read_header()) {
    ErrorDialog(XtParent(w),
                "An error occured processing nodal result data!\n"
                "No further input processing!");
    stopReading = 1;
  }
  else {
    if (i>=blocks) {
      blocks+=10;
      SETS55=(set_55*) XtRealloc( (char*) SETS55, sizeof(set_55)*blocks);
    }
    SETS55[i].header=set55.header;
    for (n=0; n<SETS55[i].header.n_int_data_val; ++n) {
      SETS55[i].header.type_spec_int_par[n] =
        set55.header.type_spec_int_par[n];
    }
    for (n=0; n<SETS55[i].header.n_real_data_val; ++n)
      SETS55[i].header.type_spec_double_par[n] =
        set55.header.type_spec_double_par[n];
    type_size=2;
    if (set55.header.data_type==2)
      type_size=1;
    SETS55[i].dat = (node_dat*)
                    XtMalloc(type_size*(N_NODES+1)*sizeof(node_dat));

    int upper=set55.header.n_data_val_per_node*type_size;
    do {
      if (!set55.read_data(&status)) {
        ErrorDialog(XtParent(w),
                    "An error occured processing nodal result data!\n"
                    "No further input processing!");
        stopReading = 1;
        i--;
      }
      else {
        for (n=0; n<upper; ++n) {
          j = nodeReverseLabels[set55.data.node_num];
          SETS55[i].dat[j].data[n] =
            (double) set55.data.node_data[n];
        }
        if ((++k % 500)==0) {
          CheckForInterrupt();
          k = 0;
        }
      }
    } while (status!=dataset::SETEmpty && !stopReading);
    sprintf(buffer, " done.\n");
    wprint(work,buffer);
    if (!set55.close_set()) {
      ErrorDialog(XtParent(w),
                  "An error occured processing nodal result data!\n"
                  "No further input processing!");
      stopReading = 1;
    }
  }
  return stopReading;
}

int ReadEResultsUnivFile(int &i, Widget w, unv *unvf, int &blocks,
                         int *nodeReverseLabels, int *elemReverseLabels) {

  dataset_56 set56;
  char       buffer[256];
  int        n,type_size,k=0,stopReading=0,status;
  int        localElementNumber;

  set56.connect(unvf);
  ++i;
  sprintf(buffer, "Reading element results no. %d:",i+1);
  wprint(work,buffer);
  if (!set56.read_header()) {
    ErrorDialog(XtParent(w),
                "An error occured processing element result data!\n"
                "No further input processing!");
    stopReading = 1;
  }
  else {
    if (i>=blocks) {
      blocks+=10;
      SETS56=(set_56*) XtRealloc( (char*) SETS56, sizeof(set_56)*blocks);
    }
    SETS56[i].header=set56.header;
    for (n=0; n<SETS56[i].header.n_int_data_val; ++n) {
      SETS56[i].header.type_spec_int_par[n] =
        set56.header.type_spec_int_par[n];
    }
    for (n=0; n<SETS56[i].header.n_real_data_val; ++n)
      SETS56[i].header.type_spec_double_par[n] =
        set56.header.type_spec_double_par[n];
    //type_size=2;
    type_size=1;
    if (set56.header.data_type==2)
      type_size=1;
    SETS56[i].dat = (elem_dat*)
                    XtCalloc(type_size*(N_ELEMENTS+1),sizeof(elem_dat));
    localElementNumber = 0;

    do {
      if (!set56.read_data(&status)) {
        if (status!=dataset::SETEmptySet) {
          ErrorDialog(XtParent(w),
                      "An errror occured reading element results!\n"
                      " Check universal file!");
          TimeoutCursors(FALSE, FALSE);
          i--;
          return 0;
        }
        else {
          i--;
          status = dataset::SETEmpty;
        }
      }
      else {
        int upper=set56.data.n_val_on_element*type_size;
        localElementNumber = 
          elemReverseLabels[set56.data.element_num];
        for (n=0; n<upper; ++n) {
          // j = set56.data.element_num;
          SETS56[i].dat[localElementNumber].data[n] =
            (double) set56.data.element_data[n];
        }
        if ((++k % 500)==0) {
          CheckForInterrupt();
          k = 0;
        }
      }
    } while (status!=dataset::SETEmpty && !stopReading);
    sprintf(buffer, " done.\n");
    wprint(work,buffer);
    if (!set56.close_set()) {
      ErrorDialog(XtParent(w),
                  "An error occured processing element result data!\n"
                  "No further input processing!");
      stopReading = 1;
    }
  }
  return stopReading;
}

// introduced by Simon Triebenbacher
int ReadUnvDimension(unv *unvf)
{
  dataset_666 set666;
  int stat = 0;

  if (!unvf) {
    return -1;
  }

  set666.connect(unvf);

  if(!set666.read_data(&stat))
    return -1;
  else
    return set666.get_dim();
}

int ReadUniversalFile(Widget w) {

  unv     unvf;
  dataset set;
  char    buffer[255];
  int     setid, maxNodeLabel,maxElemLabel;
  int     *nodeReverseLabels,*elemReverseLabels;
  int     n_blocks=10, e_blocks=10, stopReading=0, i55,i56;

  free(SETS55);
  N_SETS55 = 0;
  SETS55   = (set_55*) XtMalloc(sizeof(set_55)*n_blocks);

  free(SETS56);
  N_SETS56 = 0;
  SETS56   = (set_56*) XtMalloc(sizeof(set_56)*e_blocks);

  i55 = i56 = -1;

  TimeoutCursors(TRUE, FALSE);

  if(!unvf.open("", UNV_FILE))
    return -1;

  N_NODES    = 0;
  N_ELEMENTS = 0;
  DIM        = 0;

  if (NODES!=NULL)
    XtFree((char*) NODES);
  NODES = (point*) XtMalloc((MAX_NUM_NODES+1)*sizeof(point));
  if (ELEMENTS!=NULL)
    XtFree((char*) ELEMENTS);
  ELEMENTS = (element_data*) XtMalloc(
    (MAX_NUM_ELEMENTS+1)*sizeof(element_data));

  set.connect(&unvf);
  set.open_set(&setid);

  while (setid!=dataset::SETEOF && !stopReading) {

    switch (setid) {

      // introduced by Simon Triebenbacher
    case -666:
      DIM = ReadUnvDimension(&unvf);
      if(DIM < 0)
        return -1;
      break;

    case 781:
      maxNodeLabel =
        ReadNodesUnivFile(&nodeReverseLabels,w,&unvf);
      if (!maxNodeLabel) return 0;
      MAX_NODE_LABEL = maxNodeLabel;
      break;

    case 780:
      maxElemLabel = 
        ReadElementsUnivFile(nodeReverseLabels, w, &unvf,
                             maxNodeLabel, &elemReverseLabels);
      if (!maxElemLabel) return 0;
      break;

    case 55:
      stopReading =
        ReadNResultsUnivFile(i55,w,&unvf,n_blocks,nodeReverseLabels);
      break;

    case 56:
      stopReading =
        ReadEResultsUnivFile(i56,w,&unvf,e_blocks,
                             nodeReverseLabels,elemReverseLabels);

      break;

    default:
      break;
    }
    set.open_set(&setid);
  }


  N_SETS55=i55+1;
  N_SETS56=i56+1;

  sprintf(buffer, "Number of nodes in model: %d.\n",N_NODES);
  wprint(work,buffer);
  sprintf(buffer, "Number of elements in model: %d.\n",N_ELEMENTS);
  wprint(work,buffer);
  sprintf(buffer, "Number of nodal result sets read: %d.\n",N_SETS55);
  wprint(work,buffer);
  sprintf(buffer, "Number of element result sets read: %d.\n",N_SETS56);
  wprint(work,buffer);
  sprintf(buffer, "Setting up geometry information.\n");
  wprint(work,buffer);

  CheckForInterrupt();

  TimeoutCursors(FALSE, FALSE);

  return 0;
}
