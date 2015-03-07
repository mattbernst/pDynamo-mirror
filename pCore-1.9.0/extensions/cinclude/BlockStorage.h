/*------------------------------------------------------------------------------
! . File      : BlockStorage.h
! . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
# ifndef _BLOCKSTORAGE
# define _BLOCKSTORAGE

# include "Definitions.h"
# include "List.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The block type. */
typedef struct {
   Integer    ndata     ;
   Integer16 *indices16 ;
   Integer32 *indices32 ;
   Real      *data      ;
} Block ;

/* . The blockstorage type. */
typedef struct {
    Boolean QUNDERFLOW ;
    Integer blocksize  ;
    Integer ndata      ;
    Integer nindices16 ;
    Integer nindices32 ;
    Real    underflow  ;
    List   *blocks     ;
} BlockStorage ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern BlockStorage *BlockStorage_Allocate   ( void ) ;
extern void          BlockStorage_Data_Add   ( BlockStorage  *self, const Integer ndata, const Real *data, const Integer16 *indices16, const Integer32 *indices32 ) ;
extern void          BlockStorage_Deallocate ( BlockStorage **self ) ;
extern Block        *BlockStorage_Iterate    ( BlockStorage  *self ) ;
extern Real          BlockStorage_Size       ( BlockStorage  *self ) ;

# endif
