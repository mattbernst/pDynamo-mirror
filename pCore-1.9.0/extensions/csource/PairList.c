/*------------------------------------------------------------------------------
! . File      : PairList.c
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
!
! . Pairlists are lists of distinct pairs of integers (i,j).
!
! . Cross-pairlists are those in which one index of a pair comes from one set
! . and the other index from another set.
!
! . Self-pairlists are those in which the indices of all pairs come from the
! . same set.
!
! . Pairlists are (or should be) always sorted. The first (i) and second (j)
! . indices occur in ascending order for both cross- and self-pairlists. In
! . addition, for self-pairlists i > j.
!
! . Due to the similarity in data structures and procedures, cross- and self-
! . pairlists are handled together. Cross-pairlists are more fundamental than
! . self-pairlists however (c.f matrices and symmetric matrices).
!
!=================================================================================================================================*/

# include <stdio.h>

# include "Definitions.h"
# include "Memory.h"
# include "PairList.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer PairIndex_Compare ( const void *vterm1, const void *vterm2 ) ;

/*==================================================================================================================================
! . General pairlist procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *PairList_Allocate ( const Boolean QSELF )
{
    PairList *self = NULL ;
    self = ( PairList * ) Memory_Allocate ( sizeof ( PairList ) ) ;
    if ( self != NULL )
    {
        /* . Basic representation. */
        self->QSELF       = QSELF ;
        self->npairs      = 0 ;
        self->upperBoundi = 0 ;
        self->upperBoundj = 0 ;
        self->pairs       = List_Allocate ( ) ;
        if ( self->pairs != NULL )
        {
            self->pairs->Element_Deallocate = IndexedSelection_ListDeallocate ;
            /* . Alternative representations. */
            self->nconnectionsi   =   -1 ;
            self->connectionsi    = NULL ;
            self->connectionsj    = NULL ;
            self->numberOfRecords =   -1 ;
            self->records         = NULL ;
        }
        else PairList_Deallocate ( &self ) ;
        /* . Alternative representations. */
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear representations.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairList_ClearRepresentations ( PairList *self )
{
    if ( self != NULL )
    {
        free ( self->connectionsi ) ;
        free ( self->connectionsj ) ;
        free ( self->records      ) ;
        self->nconnectionsi   =   -1 ;
        self->connectionsi    = NULL ;
        self->connectionsj    = NULL ;
        self->numberOfRecords =   -1 ;
        self->records         = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairList_Deallocate ( PairList **self )
{
    if ( (*self) != NULL )
    {
        PairList_ClearRepresentations ( (*self) ) ;
        List_Deallocate ( &((*self)->pairs) ) ;
        free ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find an ordered pairlist from a single array containing pair indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *PairList_FromIntegerPairArray ( const Boolean QSELF, const Integer npairs, Integer *indices )
{
    PairList *self = NULL ;
    if ( ( npairs > 0 ) && ( indices != NULL ) )
    {
        auto Integer  i, iold, j, jold, m, maxint, n, p, q ;
        auto Integer *interactions ;
        /* . For self-pairlists make sure index1 is bigger than index2. */
        if ( QSELF )
        {
            for ( n = p = 0 ; p < npairs ; n += 2, p++ )
            {
                i = indices[n]   ;
                j = indices[n+1] ;
                if ( i < j ) { indices[n] = j ; indices[n+1] = i ; }
            }
        }
        /* . Sort the array. */
        qsort ( ( void * ) indices, ( size_t ) npairs, 2 * sizeof ( Integer ), ( void * ) PairIndex_Compare ) ;
        /* . Find the biggest number of interactions for an index. */
        for ( iold = -1, maxint = m = n = p = 0 ; p < npairs ; p++ )
        {
            i = indices[n] ;
            if ( i != iold )
            {
                maxint = Maximum ( maxint, m ) ;
                m = 0 ;
            }
            iold = i ;
            m++ ;
            n   += 2 ;
        }
        maxint = Maximum ( maxint, m ) ;
        /* . Allocate space. */
        interactions = Memory_Allocate_Array_Integer ( maxint ) ;
        self         = PairList_Allocate ( QSELF ) ;
        if ( ( interactions != NULL ) && ( self != NULL ) )
        {
            auto Boolean                   QOK = True       ;
            auto IndexedSelection *indexedSelection ;
            /* . Create the  pairlist without duplicates and self terms. */
            m = 0 ;
            q = 0 ;
            for ( iold = jold = -1, n = p = 0 ; p < npairs ; p++ )
            {
                i = indices[n]   ;
                j = indices[n+1] ;
                if ( i != j )
                {
                    if ( ( i != iold ) && ( m > 0 ) )
                    {
                        indexedSelection = IndexedSelection_Allocate ( iold, m, interactions ) ;
                        if ( indexedSelection == NULL ) QOK = False ;
	                List_Element_Append ( self->pairs, ( void * ) indexedSelection ) ;
                        m = 0 ;
                    }
                    if ( ( i != iold ) || ( j != jold ) )
                    {
                        interactions[m] = j ;
                        m ++ ;
                        q ++ ;
                    }
                }
                iold = i ;
                jold = j ;
                n   += 2 ;
            }
            if ( m > 0 )
            {
                indexedSelection = IndexedSelection_Allocate ( iold, m, interactions ) ;
                if ( indexedSelection == NULL ) QOK = False ;
	        List_Element_Append ( self->pairs, ( void * ) indexedSelection ) ;
            }
            /* . Finish up. */
            self->npairs      = q ;
            self->upperBoundi = indices[2 * ( npairs - 1 )] + 1 ;
            self->upperBoundj = indices[2 *   npairs - 1  ] + 1 ;
            Memory_Deallocate_Integer ( &interactions )         ;
            if ( ! QOK ) PairList_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iterate over the interactions in a list.
!---------------------------------------------------------------------------------------------------------------------------------*/
IndexedSelection *PairList_Iterate ( PairList *self )
{
    if ( self == NULL ) return NULL ;
    else                return ( IndexedSelection * ) List_Iterate ( self->pairs ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the length of the pairlist.
!---------------------------------------------------------------------------------------------------------------------------------*/
int PairList_Length ( const PairList *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->npairs ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the records representation of the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairList_MakeRecords ( PairList *self )
{
    if ( ( self != NULL ) && ( self->records == NULL ) )
    {
        auto Integer n ;
        n = PairList_NumberOfRecords ( self ) ;
        MEMORY_ALLOCATEARRAY ( self->records, n, IndexedSelection * ) ;
        if ( self->records != NULL )
        {
            auto IndexedSelection *record ;
            n = 0 ;
            List_Iterate_Initialize ( self->pairs ) ;
            while ( ( record = PairList_Iterate ( self ) ) != NULL ) { self->records[n] = record ; n += 1 ; }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of records in the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer PairList_NumberOfRecords ( PairList *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        if ( self->numberOfRecords >= 0 ) n = self->numberOfRecords ;
        else
        {
            auto IndexedSelection *record ;
            List_Iterate_Initialize ( self->pairs ) ;
            while ( ( record = PairList_Iterate ( self ) ) != NULL ) n += 1 ;
            self->numberOfRecords = n ;
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create an integer array containing the pairs.
! . The array is allocated.
!---------------------------------------------------------------------------------------------------------------------------------*/
int *PairList_ToIntegerPairArray ( PairList *self )
{
    Integer *indices = NULL ;
    if ( ( self != NULL ) && ( self->npairs > 0 ) )
    {
        indices = Memory_Allocate_Array_Integer ( 2 * self->npairs ) ;
        if ( indices != NULL )
        {
            auto Integer i, m, n = 0 ;
            auto IndexedSelection *indexedSelection ;
            /* . Loop over the pairs. */
            List_Iterate_Initialize ( self->pairs ) ;
            while ( ( indexedSelection = PairList_Iterate ( self ) ) != NULL )
            {
                i = indexedSelection->index ;
       	        for ( m = 0 ; m < indexedSelection->nindices ; m++, n += 2 )
	        {
                    indices[n]   = i ;
                    indices[n+1] = indexedSelection->indices[m] ;
                }
            }
        }
    }
    return indices ;
}

/*==================================================================================================================================
! . Cross-pairlist procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the connection representation of the pairlist.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status CrossPairList_MakeConnections ( PairList *self, const Integer upperBound )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( ! self->QSELF ) )
    {
        /* . The representation needs redoing. */
        if ( ( self->connectionsi == NULL ) || ( self->connectionsj == NULL ) || ( self->nconnectionsi < self->upperBoundi ) || ( ( upperBound >= 0 ) && ( self->nconnectionsi < upperBound ) ) )
        {
            /* . Initialization. */
            PairList_ClearRepresentations ( self ) ;
            /* . Define nconnectionsi as the bigger of the upperbound argument or of the upperbound in the list. */
            self->nconnectionsi = Maximum ( self->upperBoundi, upperBound ) ;
            /* . Allocate space. */
            self->connectionsi = Memory_Allocate_Array_Integer_Initialize ( self->nconnectionsi + 1, 0 ) ;
            self->connectionsj = Memory_Allocate_Array_Integer            ( self->npairs               ) ;
            /* . Allocation failed. */
            if ( ( self->connectionsi == NULL ) || ( self->connectionsj == NULL ) ) status = Status_OutOfMemory ;
            /* . Allocation OK. */
            else
            {
                auto Integer           i, j, m, n = 0   ;
                auto IndexedSelection *indexedSelection ;
                /* . Find the number of connections per point and fill connectionsj. */
                List_Iterate_Initialize ( self->pairs ) ;
                while ( ( indexedSelection = PairList_Iterate ( self ) ) != NULL )
                {
                    i = indexedSelection->index ;
                    self->connectionsi[i] = indexedSelection->nindices ;
       	            for ( m = 0 ; m < indexedSelection->nindices ; m++, n++ )
	            {
	                j = indexedSelection->indices[m] ;
                        self->connectionsj[n] = j ;
                    }
                }
                /* . Construct connectionsi. */
                for ( i = 0, n = 0 ; i < self->nconnectionsi ; i++ )
                {
                    m = self->connectionsi[i] ;
                    self->connectionsi[i] = n  ;
                    n += m ;
                }
                self->connectionsi[self->nconnectionsi] = n ;
                /* . Finish up. */
                status = Status_Success ;
            }
        }
        else status = Status_Success ;
    }
    return status ;
}

/*==================================================================================================================================
! . Self-pairlist procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a self-pairlist into another one.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *SelfPairList_FromSelfPairList ( PairList *self, Selection *andSelection, Selection *orSelection, const Boolean QRENUMBER )
{
    PairList *new = NULL ;
    if ( ( self != NULL ) && ( self->QSELF ) )
    {
        auto Boolean          *QAND = NULL, QNEWNUMBERS, QOK, *QOR = NULL, QORTEST ;
        auto Integer           i, *indices, j, m, n, nindices ;
        auto IndexedSelection *indexedSelection, *newIndexedSelection ;

        /* . Check the AND selection. */
        if ( andSelection == NULL )
        {
            QAND        = Memory_Allocate_Array_Boolean_Initialize ( self->upperBoundi, True ) ;
            QNEWNUMBERS = False ;
        }
        else
        {
            Selection_MakeFlags ( andSelection, self->upperBoundi ) ;
            QAND        = andSelection->flags ;
            QNEWNUMBERS = QRENUMBER           ;
            if ( QNEWNUMBERS ) Selection_MakePositions ( andSelection, self->upperBoundi ) ;
        }

        /* . Check the OR selection. */
        QORTEST = ( orSelection != NULL ) ;
        if ( QORTEST )
        {
            Selection_MakeFlags ( orSelection, self->upperBoundi ) ;
            QOR = orSelection->flags ;
        }

        /* . Create the pairlist and other temporary space. */
        new     = PairList_Allocate ( True ) ;
        indices = Memory_Allocate_Array_Integer ( self->upperBoundi - 1 ) ;

        /* . Check for a memory error. */
        QOK = ( new != NULL ) && ( indices != NULL ) && ( QAND != NULL ) ;
        if ( QOK )
        {
            /* . Iterate over the list. */
            List_Iterate_Initialize ( self->pairs ) ;
            while ( ( indexedSelection = PairList_Iterate ( self ) ) != NULL )
            {
                i = indexedSelection->index ;
                if ( QAND[i] )
                {
                    /* . Loop over the indices. */
       	            for ( m = n = 0 ; m < indexedSelection->nindices ; m++ )
	            {
                        j = indexedSelection->indices[m] ;
                        if ( QAND[j] )
                        {
                            indices[n] = j ;
                            n++ ;
                        }
                    }
                    nindices = n ;
                    /* . Apply OR test. */
                    if ( QORTEST && ( ! QOR[i] ) )
                    {
                        /* . Loop over the indices. */
       	                for ( m = n = 0 ; m < nindices ; m++ )
	                {
                            j = indices[m] ;
                            if ( QOR[j] )
                            {
                                indices[n] = j ;
                                n++ ;
                            }
                        }
                        nindices = n ;
                    }
                    /* . Renumber. */
                    if ( QNEWNUMBERS )
                    {
                        i = andSelection->positions[i] ;
                        /* . Loop over the indices. */
       	                for ( m = 0 ; m < nindices ; m++ )
	                {
                            j = indices[m] ;
                            indices[m] = andSelection->positions[j] ;
                        }
                    }
                    /* . Save the data. */
                    if ( nindices > 0 )
                    {
                        newIndexedSelection = IndexedSelection_Allocate ( i, nindices, indices ) ;
                        /* . No memory left. */
                        if ( newIndexedSelection == NULL )
                        {
                            QOK = False ;
                            break ;
                        }
                        /* . Memory left. */
                        else
                        {
  	                    List_Element_Append ( new->pairs, ( void * ) newIndexedSelection ) ;
	                    new->npairs     += nindices ;
                            new->upperBoundi = i + 1    ;
                            new->upperBoundj = indices[nindices-1] + 1 ;
                        }
                    }
                }
            }
        }
        /* . Finish up. */
        Memory_Deallocate_Integer ( &indices ) ;
        if ( andSelection == NULL ) Memory_Deallocate_Boolean ( &QAND ) ;
        if ( ! QOK ) PairList_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the connection representation of the pairlist.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status SelfPairList_MakeConnections ( PairList *self, const Integer upperBound )
{
    Status status = Status_Null ;
    if ( ( self != NULL ) && ( self->QSELF ) )
    {
        /* . The representation needs redoing. */
        if ( ( self->connectionsi == NULL ) || ( self->connectionsj == NULL ) || ( self->nconnectionsi < self->upperBoundi ) || ( ( upperBound >= 0 ) && ( self->nconnectionsi < upperBound ) ) )
        {
            auto Integer *connectionsn ;
            /* . Initialization. */
            PairList_ClearRepresentations ( self ) ;
            /* . Define nconnectionsi as the bigger of the upperbound argument or of the upperbound in the list. */
            self->nconnectionsi = Maximum ( self->upperBoundi, upperBound ) ;
            /* . Allocate space. */
            self->connectionsi = Memory_Allocate_Array_Integer_Initialize ( self->nconnectionsi + 1, 0 ) ;
            self->connectionsj = Memory_Allocate_Array_Integer            ( 2 * self->npairs           ) ;
            connectionsn       = Memory_Allocate_Array_Integer_Initialize ( self->nconnectionsi,     0 ) ;
            /* . Allocation failed. */
            if ( ( self->connectionsi == NULL ) || ( self->connectionsj == NULL ) || ( connectionsn == NULL ) )
            {
                status = Status_OutOfMemory ;
            }
            /* . Allocation OK. */
            else
            {
                auto Integer                    i, j, n          ;
                auto IndexedSelection *indexedSelection ;
                /* . Find the number of connections per point. */
                List_Iterate_Initialize ( self->pairs ) ;
                while ( ( indexedSelection = PairList_Iterate ( self ) ) != NULL )
                {
                    i  = indexedSelection->index ;
       	            for ( n = 0 ; n < indexedSelection->nindices ; n++ )
	            {
	                j = indexedSelection->indices[n] ;
                        connectionsn[i] += 1 ;
                        connectionsn[j] += 1 ;
                    }
                }
                /* . Construct connectionsi and reinitialize connectionsn. */
                for ( i = 0, n = 0 ; i < self->nconnectionsi ; i++ )
                {
                    n += connectionsn[i] ;
                    self->connectionsi[i+1] = n ;
                    connectionsn[i]         = 0 ;
                }
                /* . Fill connectionsj. */
                List_Iterate_Initialize ( self->pairs ) ;
                while ( ( indexedSelection = PairList_Iterate ( self ) ) != NULL )
                {
                    i  = indexedSelection->index ;
       	            for ( n = 0 ; n < indexedSelection->nindices ; n++ )
	            {
	                j = indexedSelection->indices[n] ;
                        self->connectionsj[self->connectionsi[i]+connectionsn[i]] = j ;
                        self->connectionsj[self->connectionsi[j]+connectionsn[j]] = i ;
                        connectionsn[i] += 1 ;
                        connectionsn[j] += 1 ;
                    }
                }
                /* . Finish up. */
                Memory_Deallocate ( connectionsn ) ;
                status = Status_Success ;
            }
        }
        else status = Status_Success ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a self-pairlist into a cross-pairlist.
! . For the moment, the connections representation is used.
! . No check is made on overlapping selections except if the flag to exclude
! . self-pairs is set.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *SelfPairList_ToCrossPairList ( PairList *self, Selection *andSelection1, Selection *andSelection2, Selection *orSelection, const Boolean QEXCLUDESELF, const Boolean QINCLUDESELF, const Boolean QRENUMBER1, const Boolean QRENUMBER2 )
{
    PairList *new = NULL ;
    if ( ( self != NULL ) && ( self->QSELF ) )
    {
        auto Boolean          *QAND1 = NULL, *QAND2 = NULL, QNEWNUMBERS1, QNEWNUMBERS2, QOK, *QOR = NULL, QORTEST ;
        auto Integer           i, index, *indices, j, m, n, nindices ;
        auto IndexedSelection *newIndexedSelection ;
        auto Status            status ;

        /* . Check the AND selections. */
        if ( andSelection1 == NULL )
        {
            QAND1        = Memory_Allocate_Array_Boolean_Initialize ( self->upperBoundi, True ) ;
            QNEWNUMBERS1 = False ;
        }
        else
        {
            Selection_MakeFlags ( andSelection1, self->upperBoundi ) ;
            QAND1        = andSelection1->flags ;
            QNEWNUMBERS1 = QRENUMBER1           ;
            if ( QNEWNUMBERS1 ) Selection_MakePositions ( andSelection1, self->upperBoundi ) ;
        }
        if ( andSelection2 == NULL )
        {
            QAND2        = Memory_Allocate_Array_Boolean_Initialize ( self->upperBoundi, True ) ;
            QNEWNUMBERS2 = False ;
        }
        else
        {
            Selection_MakeFlags ( andSelection2, self->upperBoundi ) ;
            QAND2        = andSelection2->flags ;
            QNEWNUMBERS2 = QRENUMBER2          ;
            if ( QNEWNUMBERS2 ) Selection_MakePositions ( andSelection2, self->upperBoundi ) ;
        }

        /* . Check the OR selection. */
        QORTEST = ( orSelection != NULL ) ;
        if ( QORTEST )
        {
            Selection_MakeFlags ( orSelection, self->upperBoundi ) ;
            QOR = orSelection->flags ;
        }

        /* . Create the pairlist and other temporary space. */
        new     = PairList_Allocate ( False ) ;
        status  = SelfPairList_MakeConnections  ( self, -1 ) ;
        indices = Memory_Allocate_Array_Integer ( self->upperBoundi ) ;

        /* . Check for a memory error. */
        QOK = ( new != NULL ) && ( indices != NULL ) && ( QAND1 != NULL ) && ( QAND2 != NULL ) && ( status == Status_Success ) ;
        if ( QOK )
        {
            /* . Iterate over the list. */
            for ( i = 0 ; i < Minimum ( self->nconnectionsi, self->upperBoundi ) ; i++ )
            {
                if ( QAND1[i] )
                {
                    /* . Loop over the indices. */
       	            for ( m = self->connectionsi[i], n = 0 ; m < self->connectionsi[i+1] ; m++ )
	            {
                        j = self->connectionsj[m] ;
                        if ( QAND2[j] )
                        {
                            indices[n] = j ;
                            n++ ;
                        }
                    }
                    nindices = n ;
                    /* . Apply OR test. */
                    if ( QORTEST && ( ! QOR[i] ) )
                    {
                        /* . Loop over the indices. */
       	                for ( m = n = 0 ; m < nindices ; m++ )
	                {
                            j = indices[m] ;
                            if ( QOR[j] )
                            {
                                indices[n] = j ;
                                n++ ;
                            }
                        }
                        nindices = n ;
                    }
                    /* . Exclude self-interactions. */
                    if ( QEXCLUDESELF && QAND2[i] )
                    {
                        /* . Loop over the indices. */
       	                for ( m = n = 0 ; m < nindices ; m++ )
	                {
                            j = indices[m] ;
                            if ( i != j )
                            {
                                indices[n] = j ;
                                n++ ;
                            }
                        }
                        nindices = n ;
                    }
                    /* . Include the self-interaction. */
                    if ( QINCLUDESELF && QAND2[i] )
                    {
                        indices[n] = i ;
                        nindices++ ;
                    }
                    /* . Renumber. */
                    if ( QNEWNUMBERS1 ) index = andSelection1->positions[i] ;
                    else                index = i ;
                    if ( QNEWNUMBERS2 )
                    {
                        /* . Loop over the indices. */
       	                for ( m = 0 ; m < nindices ; m++ )
	                {
                            j = indices[m] ;
                            indices[m] = andSelection2->positions[j] ;
                        }
                    }
                    /* . Save the data. */
                    if ( nindices > 0 )
                    {
                        newIndexedSelection = IndexedSelection_Allocate ( index, nindices, indices ) ;
                        if ( QINCLUDESELF ) IndexedSelection_Sort ( newIndexedSelection ) ;
                        /* . No memory left. */
                        if ( newIndexedSelection == NULL )
                        {
                            QOK = False ;
                            break ;
                        }
                        /* . Memory left. */
                        else
                        {
 	                    List_Element_Append ( new->pairs, ( void * ) newIndexedSelection ) ;
	                    new->npairs     += newIndexedSelection->nindices  ;
                            new->upperBoundi = newIndexedSelection->index + 1 ;
                            new->upperBoundj = newIndexedSelection->indices[newIndexedSelection->nindices-1] + 1 ;
                        }
                    }
                }
            }
        }
        /* . Finish up. */
        Memory_Deallocate_Integer ( &indices ) ;
        if ( andSelection1 == NULL ) Memory_Deallocate_Boolean ( &QAND1 ) ;
        if ( andSelection2 == NULL ) Memory_Deallocate_Boolean ( &QAND2 ) ;
        if ( ! QOK ) PairList_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a self-pairlist into an isolate selection container.
!---------------------------------------------------------------------------------------------------------------------------------*/
SelectionContainer *SelfPairList_ToIsolateSelectionContainer ( PairList *self, const Integer upperBound )
{
    SelectionContainer *new = NULL ;
    if ( ( self != NULL ) && ( self->QSELF ) )
    {
        auto Boolean   *QASSIGN  = NULL, QOK ;
        auto Integer   *indicesi = NULL, *indicesj = NULL, n, upper ;
        auto Status     status    ;
        auto Selection *selection ;
        /* . Get the upperbound. */
        upper = Maximum ( self->upperBoundi, upperBound ) ;
        /* . Allocate space. */
        n        = upper ;
        indicesi = Memory_Allocate_Array_Integer ( n + 1 ) ;
        indicesj = Memory_Allocate_Array_Integer ( n     ) ;
        QASSIGN  = Memory_Allocate_Array_Boolean_Initialize ( n, False ) ;
        /* . Make the connections. */
        status = SelfPairList_MakeConnections ( self, upper ) ;
/*
{
auto Integer c, i ;
printf ( "\nConnections:\n" ) ;
for ( i = 0 ; i < upper ; i++ )
{
    for ( c = self->connectionsi[i] ; c < self->connectionsi[i+1] ; c++ ) { printf ( "i, j = %6d %6d\n", i, self->connectionsj[c] ) ; }
}
}
*/
        /* . Check for memory. */
        QOK = ( indicesi != NULL ) && ( indicesj != NULL) && ( QASSIGN != NULL ) && ( status != Status_OutOfMemory ) ;
        if ( QOK )
        {
            auto Integer c, i, j, nisolates, nstart, s ;
            /* . Loop over all indices. */
            for ( n = nisolates = s = 0 ; s < upper ; s++ )
            {
                if ( ! QASSIGN[s] )
                {
                    /* . Start the new isolate. */
                    indicesj[n] = s    ;
                    nstart      = n    ;
                    QASSIGN[s]  = True ;
                    n++ ;
                    /* . Assign all indices in the new isolate. */
                    for ( i = nstart ; i < n ; i++ )
                    {
                        for ( c = self->connectionsi[indicesj[i]] ; c < self->connectionsi[indicesj[i]+1] ; c++ )
                        {
                            j = self->connectionsj[c] ;
                            if ( ! QASSIGN[j]  )
                            {
                                indicesj[n] = j    ;
                                QASSIGN[j]  = True ;
                                n++ ;
                            }
                        }
                    }
                    /* . Create the new isolate. */
                    if ( n > nstart )
                    {
                        indicesi[nisolates] = nstart ;
                        nisolates++ ;
                    }
                }
            }
            indicesi[nisolates] = n ;
            /* . Create the isolates. */
            new = SelectionContainer_Allocate ( nisolates ) ;
            if ( new != NULL )
            {
                new->QOWNER = True ;
                for ( s = 0 ; s < nisolates ; s++ )
                {
                    n = indicesi[s+1] - indicesi[s] ;
                    selection = Selection_FromIntegerArray ( n, &(indicesj[indicesi[s]]) ) ;
                    Selection_Sort ( selection ) ;
                    new->items[s]   = selection ;
                    new->upperBound = Maximum ( Selection_UpperBound ( selection ), new->upperBound ) ;
                }
            }
            else QOK = False ;
        }
        /* . Finish up. */
        Memory_Deallocate_Integer ( &indicesi ) ;
        Memory_Deallocate_Integer ( &indicesj ) ;
        Memory_Deallocate_Boolean ( &QASSIGN  ) ;
        if ( ! QOK ) SelectionContainer_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the upper bound.
! . For self-pairlists i > j.
!---------------------------------------------------------------------------------------------------------------------------------*/
int SelfPairList_UpperBound ( PairList *self )
{
    if ( ( self == NULL ) || ( ! self->QSELF ) ) return 0 ;
    else                                         return self->upperBoundi ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
static Integer PairIndex_Compare ( const void *vterm1, const void *vterm2 )
{
    Integer   i ;
    Integer  *term1, *term2 ;
    term1 = ( Integer * ) vterm1 ;
    term2 = ( Integer * ) vterm2 ;
         if ( term1[0] < term2[0] ) i = -1 ;
    else if ( term1[0] > term2[0] ) i =  1 ;
    else if ( term1[1] < term2[1] ) i = -1 ;
    else if ( term1[1] > term2[1] ) i =  1 ;
    else i = 0 ;
    return i ;
}
