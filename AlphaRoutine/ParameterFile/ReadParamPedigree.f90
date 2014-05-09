subroutine ReadParamPedigree(LFileNameC)

use Pedigree

implicit none

integer :: i
integer :: LEndPedLineI,LLineWithPedigreeFileKeywordC,LStartPedLineI,LSearchStartI,LSeachEndI,LSubModeI,LNumOfColumnI
character(len=2000) :: LCurrentWordSearchC,LDumC,LFileNameC,LKeyWordC
logical :: LEndPedPresentL,LPedigreeFileKeywordPresentL,LStartPedPresentL

!Start identify line number of '@StartPedigree@' keyword
	LSubModeI=1
	LKeyWordC='@StartPedigree@'
	call CheckKeyWordPresence(LKeyWordC,LFileNameC,LSearchStartI,LSeachEndI,LStartPedLineI,LStartPedPresentL,LSubModeI)
!End Identify line number of '@StartPedigree@' keyword

!Start identify line number of 'EndPedigree@' keyword
	LSubModeI=1
	LKeyWordC='@EndPedigree@'
	call CheckKeyWordPresence(LKeyWordC,LFileNameC,LSearchStartI,LSeachEndI,LEndPedLineI,LEndPedPresentL,LSubModeI)
!End identify line number of 'EndPedigree@' keyword

!Start check validity of '@StartPedigree@' and '@EndPedigree@'
	if ((LStartPedPresentL==.true.).and.(LEndPedPresentL==.true.)) then
	
	else
		write(*,*) 'Error in the validity of @StartPedigree@ and @EndPedigree@'
		write(*,*) 'Check pedigree commands in spec file'
		stop
	endif
!End check validity of '@StartPedigree@' and '@EndPedigree@'	
	
!Start PedigreeFile keyword and file name search
	LSearchStartI=LStartPedLineI+1
	LSeachEndI=LEndPedLineI-1
	LKeyWordC='FileName'
	LSubModeI=2
	call CheckKeyWordPresence(LKeyWordC,LFileNameC,LSearchStartI,LSeachEndI,LLineWithPedigreeFileKeywordC,LPedigreeFileKeywordPresentL,LSubModeI)

	
	!Start check validity of PedigreeFile keyword
		if (LPedigreeFileKeywordPresentL==.true.) then
	
		else
			write(*,*) 'Error in the validity of PedigreeFile keyword'
			write(*,*) 'Check pedigree commands in spec file'
			stop
		endif
	!End check validity of PedigreeFile keyword	

	!Start read pedigree file name
	GPedFileNameC='ll'
	print*, GPedFileNameC
		!LNumOfColumnI=2
		!call GetRequiredCommandFileInfo(LFileNameC,LLineWithPedigreeFileKeywordC,LNumOfColumnI,GPedFileNameC)
		!print*, trim(GPedFileNameC)

	!End read pedigree file name
	
	close(1)
!End PedigreeFile keyword and file name search

end subroutine ReadParamPedigree

