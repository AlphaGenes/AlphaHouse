program test

  use HashModule
  implicit none

  type(DictStructure) :: t
  character(len=100000)  :: f
  integer :: i

    t = DictStructure()
    do i = 1, 100000
    write(f, '(a2,i)') "hi",i
      call t%addKey(f, i)
    enddo
    call t%destroy()


    t = DictStructure()
    do i = 1, 100000
    write(f, '(a2,i)') "hi",i
      call t%addKey(f, i)
    enddo
    call t%destroy()
  end program test
