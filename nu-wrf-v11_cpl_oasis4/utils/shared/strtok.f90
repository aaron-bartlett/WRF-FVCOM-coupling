   CHARACTER*255 FUNCTION strtok (source_string, delimiters)

! Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!         and the delimiter list used to tokenize SOURCE_STRING in DELIMITERS.
!         then, if the returned value is not equal to CHAR(0), keep calling until
!         it is with SOURCE_STRING set to CHAR(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning CHAR(0). 
!
! Input:  source_string =   Source string to tokenize. 
!         delimiters    =   delimiter string.
!                      Used to determine the beginning/end of each token in a string.
!
! Output: strtok()
!
! LIMITATIONS:
!   can not be called with a different string until current string is totally processed, even from different procedures
!   input string length limited to set size
!   function returns fixed 255 character length
!   length of returned string not given

      CHARACTER(len=*),intent(in)  :: source_string
      CHARACTER(len=*),intent(in)  :: delimiters

      CHARACTER(len=255),save :: saved_string
      INTEGER,save :: isaved_start  ! points to beginning of unprocessed data
      INTEGER,save :: isource_len   ! length of original input string

      INTEGER :: ibegin        ! beginning of token to return
      INTEGER :: ifinish       ! end of token to return

! initialize stored copy of input string and pointer into input string on first call
      IF (source_string(1:1) .NE. CHAR(0)) THEN
         isaved_start = 1                 ! beginning of unprocessed data
         saved_string = source_string     ! save input string from first call in series
         isource_len = LEN(saved_string)  ! length of input string from first call
      ENDIF

      ibegin = isaved_start

      DO
         IF ( (ibegin .LE. isource_len) .AND. (INDEX(delimiters,saved_string(ibegin:ibegin)) .NE. 0)) THEN
            ibegin = ibegin + 1
         ELSE
            EXIT
         ENDIF
      ENDDO

      IF (ibegin .GT. isource_len) THEN
         strtok = CHAR(0)
         RETURN
      ENDIF

      ifinish = ibegin

      DO
         IF ((ifinish .LE. isource_len) .AND.  (INDEX(delimiters,saved_string(ifinish:ifinish)) .EQ. 0)) THEN
            ifinish = ifinish + 1
         ELSE
            EXIT
         ENDIF
      ENDDO

      !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
      strtok = saved_string(ibegin:ifinish-1)
      isaved_start = ifinish

   END FUNCTION strtok
