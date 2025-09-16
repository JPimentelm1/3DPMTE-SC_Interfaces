C --- This is the start of the include file ---
C
C
C     Declare the type and size of every variable in the block
      REAL*8 ALLSEv, ALLWKv
      REAL*8 toler
C     Define a named common block for shared variables
      COMMON /MySharedVars/ ALLSEv, toler, ALLWKv
      INTEGER, PARAMETER :: MAX_NODES = 10000
	  INTEGER :: recordedNodes(MAX_NODES), numRecorded = 0
      LOGICAL :: nodeRecorded(MAX_NODES) = .FALSE.
C		
C --- END