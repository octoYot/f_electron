module precision
USE f_e_data   ! GLOBAL  io
! Precision:
! All real variables defaulted to double precision
   integer, parameter  :: realkind = selected_real_kind(p=13,r=200)  
end module precision


module f_e_interpreter
!
! This module interprets the function, builds it to be evaluated next
!
   use precision 
   implicit none

   character(len=10),  dimension(:), allocatable, private :: varnames 
   character(len=255), dimension(:), allocatable, private :: stokens
   integer,            dimension(:), allocatable, private :: operations
   integer,                                       private :: n 
   integer,                                       private :: ntokens = 0 
   character,          dimension(:), pointer,  private :: opaddsub   !Operator
   integer,                                    private :: isaddsub = 1
   character,          dimension(:), pointer,  private :: opmuldiv   !Operator
   integer,                                    private :: ismuldiv = 1
   character(len=255),                         private :: toke
   integer,                                    private :: itoke = 1
   integer,                                    private :: ioperations = 1
   integer,                                    private :: numberk = 1
   real(kind=realkind), dimension(:), pointer,  private :: pdata
   real(kind=realkind), dimension(:), pointer,  private :: number
   character(len=5),           public :: statusflagparser = 'ok'

contains

!
!------------------------------------------------------------------------------
!
     subroutine init(func, variablenames, statusflag)
!
!  This subroutine shifts all characters of the function 
!  expression to lowercase and converts exponents signals ** to ^
!
      character(len=*),                intent(inout)  :: func
      character(len=10), dimension(:), intent(inout)  :: variablenames 
      character(len=26)         :: lower = 'abcdefghijklmnopqrstuvwxyz'
      character(len=26)         :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      integer            :: i, k, funclen
      character(len=5)      , intent(out)  :: statusflag

      funclen = len_trim(func)
      if (funclen == 0) then
        write(io,'(/,"****FATAL: There is no function defined for a parameter/constant",  &
                   /,"Make sure there is a space as a delimiter after the parameter/constant name")') 
        STOP
      end if
!  detects errors
      call identify(func)
      call convert_b(func)
   
!  Shift all characters to lowercase and converts ** to ^
      do i = 1, funclen
         k = index(upper,func(i:i)) 
         if ( k /= 0) then
            func(i:i) = lower(k:k)
         end if
         k = index(func,'**')
         if (k /= 0) then
           func = func(:k-1) // '^' // func(k+2:)
         end if
         k = index(func,'=')
         if (k /= 0) then
           write(io,'(/,"****FATAL: You cannot use ""="" in the definition of a parameter or constant")') 
           STOP
         end if
      end do
      call blanks(func)
   
      call recog_variables (func, variablenames)
      statusflag = statusflagparser

      end subroutine init
!
!------------------------------------------------------------------------------
!
      subroutine recog_variables (func, variablenames)
!
! This subroutine recognizes the variables and set their values
!
      character(len=10), dimension(:), intent(in)    :: variablenames 
      character(len=*),                intent(inout) :: func

      n = size(variablenames)
      allocate(varnames(n))
      varnames = variablenames
      call tokens_analyzer (func)

      end subroutine recog_variables
!
!------------------------------------------------------------------------------
!
      subroutine tokens_analyzer (func)
!
! This subroutine scans the func string storing its basic elements
!
      character(len=*), intent(in) :: func
      character(len=11)            :: numbers = '.0123456789'
      character(len=26)            :: chars = 'abcdefghijklmnopqrstuvwxyz'
      character(len=7)             :: operators = '+-*/^()'
      integer           :: k = 1, i = 1 
      integer           :: irightbrackets = 1, ileftbrackets = 1 
      logical           :: status = .true. 

      k = 1
      i = 1
      irightbrackets = 1
      ileftbrackets = 1 
      ntokens = 0 
      status = .true.

      do while (k <= len_trim(func))
! It's a variable, or function name
        if (index(chars,func(k:k)) /= 0) then
          status = .true.
          do while (status)
            if (index(operators, func(k+1:k+1)) == 0 .and. func(k+1:k+1) /= ' ') then
              k = k + 1
            else
              status = .false.
              k = k + 1
            end if
          end do
          ntokens = ntokens + 1
      
!  It's a number
        else if (index(numbers,func(k:k)) /= 0) then
          status = .true.
          do while (status)
            if ((index(operators, func(k+1:k+1)) == 0 .and. func(k+1:k+1) /= ' ') .or. func(k+1:k+1) == 'e' .or. func(k+1:k+1) == 'd') then
              k = k + 1
            else
              if(func(k:k) == 'e' .or. func(k:k) == 'd') then
                k = k + 1
              else
                status = .false.
                k = k + 1
              end if
            end if
          end do
          ntokens = ntokens + 1
         
!  It's an operator or delimitator
        else 
          k = k + 1
          ntokens = ntokens + 1
        end if  
      end do ! do while (k <= len_trim(func))

      allocate(stokens(ntokens))
      k = 1
      i = 1

      do while (k <= len_trim(func))
!  It's a variable, or function name
        if (index(chars,func(k:k)) /= 0) then
          stokens(i) = func(k:k)
          status = .true.
          do while (status)
            if (index(operators, func(k+1:k+1)) == 0 .and. func(k+1:k+1) /= ' ') then
              stokens(i) = trim(stokens(i)) // func(k+1:k+1)
              k = k + 1
            else
              status = .false.
              k = k + 1
              i = i + 1
            end if
          end do
            
!  It's a number
        else if (index(numbers,func(k:k)) /= 0) then
          stokens(i) = func(k:k)
          status = .true.
          do while (status)
            if ((index(operators, func(k+1:k+1)) == 0 .and. func(k+1:k+1) /= ' ') .or. func(k+1:k+1) == 'e' .or. func(k+1:k+1) == 'd') then
              stokens(i) = trim(stokens(i)) // func(k+1:k+1)
              k = k + 1
            else
              if(func(k:k) == 'e' .or. func(k:k) == 'd') then
                stokens(i) = trim(stokens(i)) // func(k+1:k+1)
                k = k + 1
              else
                status = .false.
                i = i + 1
                k = k + 1
              end if
            end if
          end do
         
!  It's an operator or delimitator
        else 
          stokens(i) = func(k:k)
          if(stokens(i) == '(')then
            irightbrackets = irightbrackets + 1
          else if(stokens(i) == ')') then
            ileftbrackets = ileftbrackets + 1
          end if   
          i = i + 1
          k = k + 1
        end if  
      end do  ! do while (k <= len_trim(func))

      if (irightbrackets /= ileftbrackets) then
        statusflagparser = 'error'
        return
      end if  
      
      itoke = 1
      isaddsub = 1
      ismuldiv = 1
      ioperations = 1
      numberk = 1
      toke = stokens(itoke)
      allocate(opaddsub(2))
      allocate(opmuldiv(2))
      allocate(number(ntokens))
      allocate(pdata(ntokens))
      allocate(operations(ntokens))

      call add_sub()
      ioperations = ioperations - 1

      end subroutine tokens_analyzer
!
!------------------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  The following subroutines call themselves recursively
!  to build the expression to be parsed based on an algorithm
!  called Recursive Descendent Parsing
!**** must compile with -frecursive option **** MJR ****
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
     subroutine add_sub ()
!
! Called by subroutine tokens_analyzer(), brackets()
! 
      call mul_div ()

      do while (trim(toke) == '+' .or. trim(toke) == '-') 
        opaddsub(isaddsub) = trim(toke)
        isaddsub = isaddsub + 1
        itoke = itoke + 1
        toke = stokens(itoke)
        call mul_div()
         
        selectcase(opaddsub(isaddsub-1))
        case('+')
          isaddsub = isaddsub - 1
          operations(ioperations) = 3
          ioperations = ioperations + 1
      
        case('-')
          isaddsub = isaddsub - 1
          operations(ioperations) = 4
          ioperations = ioperations + 1
        end select
      end do

      end subroutine add_sub
!
!------------------------------------------------------------------------------
!
      subroutine mul_div ()
!
! Called by subroutine add_sub()
!
      call unary()

      do while (trim(toke) == '*' .or. trim(toke) == '/') 
        opmuldiv(ismuldiv) = trim(toke)
        ismuldiv = ismuldiv + 1
        itoke = itoke + 1
        toke = stokens(itoke)
        call unary()

        selectcase(opmuldiv(ismuldiv-1))
        case('*')
          ismuldiv = ismuldiv - 1 
          operations(ioperations) = 5
          ioperations = ioperations + 1   
        case('/')
          ismuldiv = ismuldiv - 1
          operations(ioperations) = 6
          ioperations = ioperations + 1
        end select
      end do

      end subroutine mul_div
!
!------------------------------------------------------------------------------
!
      subroutine unary()
!
! Called by subroutine mul_div()
!   
     if (trim(toke) == '-') then
       itoke = itoke + 1
       toke = stokens(itoke)
       call pow()
       operations(ioperations) = 2
       ioperations = ioperations + 1
     else if (trim(toke) == '+') then
       itoke = itoke + 1
       toke = stokens(itoke)
       call pow()
     else
       call pow()
     end if

     end subroutine unary
!
!------------------------------------------------------------------------------
!
      subroutine pow ()
!
! Called by subroutine unary()
!
      call functions()

      if (trim(toke) == '^') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call functions()
        operations(ioperations) = 7
        ioperations = ioperations + 1
      end if

      end subroutine pow
!
!------------------------------------------------------------------------------
!
     subroutine functions ()
!
! Called by subroutine pow()
! 
      if (trim(toke) == 'sin') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 8
        ioperations = ioperations + 1
      
      else if(trim(toke) == 'cos') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 9
        ioperations = ioperations + 1
   
      else if(trim(toke) == 'tan') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 10
        ioperations = ioperations + 1
      
      else if(trim(toke) == 'asin') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 11
        ioperations = ioperations + 1
      
      else if(trim(toke) == 'acos') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 12
        ioperations = ioperations + 1
      
      else if(trim(toke) == 'atan') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 13
        ioperations = ioperations + 1
      
      else if(trim(toke) == 'sinh') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 14
        ioperations = ioperations + 1
      
      else if(trim(toke) == 'cosh') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()  
        operations(ioperations) = 15
        ioperations = ioperations + 1
         
      else if(trim(toke) == 'tanh') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 16
        ioperations = ioperations + 1
         
      else if(trim(toke) == 'sind') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 17
        ioperations = ioperations + 1

      else if(trim(toke) == 'cosd') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 18
        ioperations = ioperations + 1
         
      else if(trim(toke) == 'tand') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()  
        operations(ioperations) = 19
        ioperations = ioperations + 1
         
      else if (trim(toke) == 'log') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 20
        ioperations = ioperations + 1
         
      else if (trim(toke) == 'log10') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 21
        ioperations = ioperations + 1
      
      else if (trim(toke) == 'nint') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 22
        ioperations = ioperations + 1
         
      else if (trim(toke) == 'anint') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 23
        ioperations = ioperations + 1
        
      else if (trim(toke) == 'aint') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 24
        ioperations = ioperations + 1
         
      else if (trim(toke) == 'exp') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 25
        ioperations = ioperations + 1
      
      else if (trim(toke) == 'sqrt') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 26
        ioperations = ioperations + 1
         
      else if (trim(toke) == 'abs') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 27
        ioperations = ioperations + 1
        
      else if (trim(toke) == 'floor') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call brackets()
        operations(ioperations) = 28
        ioperations = ioperations + 1
        
      else
        call brackets()
      end if
   
      end subroutine functions
!
!------------------------------------------------------------------------------
!
      subroutine brackets()
!
! Called by subroutine functions()
!
      if (trim(toke) == '(') then
        itoke = itoke + 1
        toke = stokens(itoke)
        call add_sub()
        if (trim(toke) /= ')') then
          statusflagparser =  'error'
          return
        end if     
        if (itoke < ntokens) then
          itoke = itoke + 1
          toke = stokens(itoke)
        end if
        if (trim(toke) == '(') then
          statusflagparser =  'error'
          return
        end if 
      else
        call recog_vars ()
      end if
   
      end subroutine brackets
!
!------------------------------------------------------------------------------
!
      subroutine recog_vars()
!
!  Called from subroutine brackets.
!
      integer          :: i
      integer          :: ierror
      character(len=7) :: operators = '+-*/^()'
   
!  Expression has an error
      if (index(operators, trim(toke)) /= 0) then
       statusflagparser = 'error'
       return 
      end if 
!   
      do i = 1, n
!  It's a variable
        if (trim(toke) == varnames(i)) then
          operations(ioperations) = 28+i
          ioperations = ioperations + 1
          if (itoke < ntokens) then
            itoke = itoke + 1
            toke = stokens(itoke)
          end if
          return
        end if
      end do 

!  It's a number
      toke = trim(toke)
      read(toke, *, iostat = ierror) number(numberk)
      if (ierror /= 0) then
        statusflagparser = 'error'
        return
      else
        operations(ioperations) = 1
        ioperations = ioperations + 1
        if (itoke < ntokens) then
          itoke = itoke + 1
          toke = stokens(itoke)
        end if
        numberk = numberk + 1
      end if
!  
      end subroutine recog_vars
!
!------------------------------------------------------------------------------
!
      function evaluate (vars) result (answer)
!
! This function will evaluate the expression with the variable values supplied.
!
      real(kind = realkind), dimension(:), intent(in)  :: vars
      real(kind = realkind)        :: answer
      integer            :: st = 0
      integer            :: dt = 1
      integer            :: i  

      st = 0
      dt = 1

      do i = 1, ioperations 
        select case(operations(i))
        case (1)
          st = st + 1
          pdata(st) = number(dt)
          dt = dt + 1
        case (2)
          pdata(st) = - pdata(st)
        case (3)
          pdata(st-1) = pdata(st-1) + pdata(st)
          st = st - 1 
        case (4)
          pdata(st-1) = pdata(st-1) - pdata(st)
          st = st - 1 
        case (5)
          pdata(st-1) = pdata(st-1) * pdata(st)
          st = st - 1 
        case (6)
          pdata(st-1) = pdata(st-1) / pdata(st)
          st = st - 1 
        case (7)
          pdata(st-1) = pdata(st-1) ** pdata(st)
          st = st - 1
        case (8)
          pdata(st) = sin(pdata(st))
        case (9)
          pdata(st) = cos(pdata(st))
        case (10)
          pdata(st) = tan(pdata(st))
        case (11)
          pdata(st) = asin(pdata(st))
        case (12)
          pdata(st) = acos(pdata(st))
        case (13)
          pdata(st) = atan(pdata(st))
        case (14)
          pdata(st) = sinh(pdata(st))
        case (15)
          pdata(st) = cosh(pdata(st))
        case (16)
          pdata(st) = tanh(pdata(st))
        case (17)
          pdata(st) = sin(pdata(st)*atan(1.0d+00)/45.0d+0)    ! =sind(pdata(st))
        case (18)
          pdata(st) = cos(pdata(st)*atan(1.0d+00)/45.0d+0)    ! =cosd(pdata(st))
        case (19)
          pdata(st) = tan(pdata(st)*atan(1.0d+00)/45.0d+0)    ! =tand(pdata(st))
        case (20)
         pdata(st) = log(pdata(st))
        case (21)
          pdata(st) = log10(pdata(st))
        case (22)
          pdata(st) = nint(pdata(st))
        case (23)
          pdata(st) = anint(pdata(st))
        case (24)
          pdata(st) = aint(pdata(st))
        case (25)
          pdata(st) = exp(pdata(st))
        case (26)
          pdata(st) = sqrt(pdata(st))
        case (27)
          pdata(st) = abs(pdata(st))
        case (28)
          pdata(st) = floor(pdata(st))
        case default
          st = st + 1
          pdata(st) = vars(operations(i)-28)
        end select
      end do
! 
      answer = pdata(1)
!
      end function evaluate
!
!------------------------------------------------------------------------------
!
      function evaluate_details (vars) result (answer)
!
!  This function will evaluate the expression supplied
!  Alternative to the FUNCTION evaluate() above, not used.
!
      real(kind = realkind), dimension(:), intent(in)  :: vars
      real(kind = realkind)        :: answer
      integer            :: st = 0
      integer            :: dt = 1
      integer            :: i  
!
      st = 0
      dt = 1

      do i = 1, ioperations 
        select case(operations(i))
         
        case (1)
          st = st + 1
          pdata(st) = number(dt)
          dt = dt + 1
         
        case (2)
          pdata(st) = - pdata(st)
         
        case (3)
          pdata(st-1) = pdata(st-1) + pdata(st)
          st = st - 1 
         
        case (4)
          pdata(st-1) = pdata(st-1) - pdata(st)
          st = st - 1 
         
        case (5)
          pdata(st-1) = pdata(st-1) * pdata(st)
          st = st - 1 
         
        case (6)
          if(abs(pdata(st)) < 1.0e-30) then
            answer = -7.093987e-35
            return
          end if
          pdata(st-1) = pdata(st-1) / pdata(st)
          st = st - 1 
         
        case (7)
          if(pdata(st-1) < 0.0 .AND. (pdata(st)-int(pdata(st))) /= 0.0) then
            answer = -7.093987e-35
            return
          end if
          if(pdata(st)*log(abs(pdata(st-1)+1.0E-15)) > 65.0) then
            answer = -7.093987e-35
            return
          end if
          pdata(st-1) = pdata(st-1) ** pdata(st)
          st = st - 1
         
        case (8)
          pdata(st) = sin(pdata(st))
         
        case (9)
          pdata(st) = cos(pdata(st))
         
        case (10)
          if((abs(pdata(st)) > 89.99*3.141593/180. .and. abs(pdata(st)) < 90.01*3.141593/180)&
          .or. (abs(pdata(st)) > 269.99*3.141593/180. .and. abs(pdata(st)) < 270.01*3.141593/180)) then
             answer = -7.093987e-35
             return
          end if
          pdata(st) = tan(pdata(st))
         
        case (11)
          if(abs(pdata(st)) > 1.0) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = asin(pdata(st))
         
        case (12)
          if(abs(pdata(st)) > 1.0) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = acos(pdata(st))
         
        case (13)
          if(abs(pdata(st)) > 1.0e+10) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = atan(pdata(st))
         
        case (14)
          if(pdata(st) > 60) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = sinh(pdata(st))
         
        case (15)
          if(pdata(st) > 60) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = cosh(pdata(st))
      
        case (16)
          pdata(st) = tanh(pdata(st))
      
        case (17)
          pdata(st) = sin(pdata(st)*atan(1.0d+00)/45.0d+0)  ! =sind(pdata(st))
      
        case (18)
          pdata(st) = cos(pdata(st)*atan(1.0d+00)/45.0d+0)  ! =cosd(pdata(st))
      
        case (19)
          pdata(st) = tan(pdata(st)*atan(1.0d+00)/45.0d+0)  ! =tand(pdata(st))
      
        case (20)
          if(pdata(st) <= 1.0e-15) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = log(pdata(st))
      
        case (21)
          if(pdata(st) <= 1.0e-15) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = log10(pdata(st))
      
        case (22)
          pdata(st) = nint(pdata(st))
      
        case (23)
          pdata(st) = anint(pdata(st))
      
        case (24)
          pdata(st) = aint(pdata(st))
      
        case (25)
          if(pdata(st) > 55) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = exp(pdata(st))
      
        case (26)
          if(pdata(st) < 0) then
            answer = -7.093987e-35
            return
          end if
          pdata(st) = sqrt(pdata(st))
      
        case (27)
          pdata(st) = abs(pdata(st))
      
        case (28)
          pdata(st) = floor(pdata(st))
      
        case default
          st = st + 1
          pdata(st) = vars(operations(i)-28)
      
        end select
      
        if(abs(pdata(st)) > 1.0d+60) then
          answer = -7.093987e-35
          return
        end if
   
      end do
   
      answer = pdata(1)
!
      end function evaluate_details
!
!------------------------------------------------------------------------------
!
      subroutine destroyfunc()
      
      if (allocated(stokens)) then
        deallocate(stokens)
      end if
      if (associated(opaddsub)) then
        deallocate(opaddsub)
      end if
      if (associated(opmuldiv)) then
        deallocate(opmuldiv)
      end if
      if (associated(number)) then
        deallocate(number)
      end if
      if (associated(pdata)) then
        deallocate(pdata)
      end if
      if (allocated(operations)) then
        deallocate(operations)
      end if
      if (allocated(varnames)) then
        deallocate(varnames)
      end if
      statusflagparser = 'ok'
 
      end subroutine destroyfunc
!
!------------------------------------------------------------------------------
!
      recursive subroutine blanks(func)
!
!  This subroutine removes unnecessary blank spaces.
!  Called from subroutine init()
!
      character(len=*), intent(inout)  :: func
      integer        :: k
   
      func = adjustl(func)
      k = index(trim(func), ' ')
      if (k /= 0) then
        func = func(:k-1) // func(k+1:)
        call blanks(func)
      end if
   
   end subroutine blanks
!
!------------------------------------------------------------------------------
!MJR module end move from here
!
      subroutine identify(func)
      character(len=*), intent(inout)  :: func
!MJR      character (255) func
      character (36) variav
!
      integer nchar,i,j
!

      variav = '0123456789abcdefghijklmnopqrstuvwxyz'

      nchar = len(trim(func))

      do i = 1,nchar
        if(func(i:i) == 'A')func(i:i) = 'a'
        if(func(i:i) == 'B')func(i:i) = 'b'
        if(func(i:i) == 'C')func(i:i) = 'c'
        if(func(i:i) == 'D')func(i:i) = 'd'
        if(func(i:i) == 'E')func(i:i) = 'e'
        if(func(i:i) == 'F')func(i:i) = 'f'
        if(func(i:i) == 'G')func(i:i) = 'g'
        if(func(i:i) == 'H')func(i:i) = 'h'
        if(func(i:i) == 'I')func(i:i) = 'i'
        if(func(i:i) == 'J')func(i:i) = 'j'
        if(func(i:i) == 'K')func(i:i) = 'k'
        if(func(i:i) == 'L')func(i:i) = 'l'
        if(func(i:i) == 'M')func(i:i) = 'm'
        if(func(i:i) == 'N')func(i:i) = 'n'
        if(func(i:i) == 'O')func(i:i) = 'o'
        if(func(i:i) == 'P')func(i:i) = 'p'
        if(func(i:i) == 'Q')func(i:i) = 'q'
        if(func(i:i) == 'R')func(i:i) = 'r'
        if(func(i:i) == 'S')func(i:i) = 's'
        if(func(i:i) == 'T')func(i:i) = 't'
        if(func(i:i) == 'U')func(i:i) = 'u'
        if(func(i:i) == 'V')func(i:i) = 'v'
        if(func(i:i) == 'W')func(i:i) = 'w'
        if(func(i:i) == 'X')func(i:i) = 'x'
        if(func(i:i) == 'Y')func(i:i) = 'y'
        if(func(i:i) == 'Z')func(i:i) = 'z'
        if(func(i:i) == ',') then  !  MJR: don't accept commas
          func = 'erro'
          return
        end if
      end do

      if(func(nchar:nchar) == '-' .or. func(nchar:nchar) == '+' .or. func(nchar:nchar) == '/' .or. func(nchar:nchar) == '*') then
        func = 'erro'
        return
      end if
    
      if(func(1:1) == '*' .or. func(1:1) == '/') then
        func = 'erro'
        return
      end if
    
      do i = 1, nchar-1
        if(func(i:i+1) == '--' .or. func(i:i+1) == '-+' .or. func(i:i+1) == '-/' .or. func(i:i+1) == '-*') func = 'erro'
        if(func(i:i+1) == '+-' .or. func(i:i+1) == '++' .or. func(i:i+1) == '+/' .or. func(i:i+1) == '+*') func = 'erro'
        if(func(i:i+1) == '*-' .or. func(i:i+1) == '*+' .or. func(i:i+1) == '*/') func = 'erro'
        if(func(i:i+1) == '/-' .or. func(i:i+1) == '/+' .or. func(i:i+1) == '//' .or. func(i:i+1) == '/*') func = 'erro'
      end do
      if(trim(func) == 'erro') return
      
      do i = 1, nchar-1
        do j = 1, 36
          if(func(i:i+1) == ')'//variav(j:j)) func = 'erro'
        end do
      end do
      if(trim(func) == 'erro') return

      do i = 1, nchar-1
        do j = 1, 36
          if(func(i:i) == '0' .or. func(i:i) == 'n' .or. func(i:i) == 's' .or. func(i:i) == 'h' .or. func(i:i) == 'd' .or. func(i:i) == 'g' .or. func(i:i) == 't' .or. func(i:i) == 'p' .or. func(i:i) == 'r') then
            !   valid funtions always end in the above letters.
          else
            if(func(i:i+1) == variav(j:j)//'(') func = 'erro'
          end if
        end do
      end do
      if(trim(func) == 'erro') return
      
      
      if(nchar >= 5) then
        do i = 1,nchar-4

          if(func(i:i+4) == 'log10') then
            j = i+5
  19        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 19
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+4) == 'anint') then
            j = i+5
  20        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 20
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+4) == 'floor') then
            j = i+5
  21        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 21
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

        end do  ! do i = 1,nchar-4
      end if  ! if(nchar >= 5)

      if(nchar >= 4) then
        do i = 1,nchar-3

          if(func(i:i+3) == 'asin') then
            j = i+4
  7         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 7
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'acos') then
            j = i+4
  8         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 8
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'atan') then
            j = i+4
  9         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 9
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'sinh') then
            j = i+4
  10        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 10
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'cosh') then
            j = i+4
  11        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 11
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'tanh') then
            j = i+4
  12        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 12
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'sind') then
            j = i+4
  13        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 13
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'cosd') then
            j = i+4
  14        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 14
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'tand') then
            j = i+4
  15        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 15
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'nint') then
            j = i+4
  16        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 16
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'aint') then
            j = i+4
  17        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 17
            end if 
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+3) == 'sqrt') then
            j = i+4
  18        continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 18
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

        end do  ! do i = 1,nchar-3
      end if  ! if(nchar >= 4) then


      if(nchar >= 3) then
        do i = 1,nchar-2

          if(func(i:i+2) == 'sin') then
            j = i+3
  1         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 1
            end if
            if(func(j:j) == 'd' .or. func(j:j) == 'h') goto 51
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
  51        continue
          end if

          if(func(i:i+2) == 'cos') then
            j = i+3
  2         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 2
            end if
            if(func(j:j) == 'd' .or. func(j:j) == 'h') goto 52
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
  52        continue
          end if

          if(func(i:i+2) == 'tan') then
            j = i+3
  3         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 3
            end if
            if(func(j:j) == 'd' .or. func(j:j) == 'h') goto 53
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
  53        continue
          end if

          if(func(i:i+2) == 'log') then
            j = i+3
  4         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 4
            end if
            if(j < (nchar-1) .and. func(j:j+1) == '10') goto 54
            if(func(j:j) /= '(') then
               func = 'erro'
              return
            end if
  54        continue
          end if

          if(func(i:i+2) == 'exp') then
            j = i+3
  5         continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 5
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if

          if(func(i:i+2) == 'abs') then
            j = i+3
 6          continue
            if(j >= nchar) then
              func = 'erro'
              return
            end if
            if(func(j:j) == ' ') then
              j = j + 1
              goto 6
            end if
            if(func(j:j) /= '(') then
              func = 'erro'
              return
            end if
          end if 

        end do ! do i = 1,nchar-2
      endif  ! if(nchar >= 3) then

      return
      end subroutine identify
!
!------------------------------------------------------------------------------
!
      subroutine convert_b(text)
      character(len=*), intent(inout)  :: text
!MJR      character (255) text
!MJR      INTENT(INOUT)::text
      integer ilength,k,item

      ilength = len(trim(text))

      do k = 1,ilength
        if(text(k:k) == '[' .or. text(k:k) == '{') text(k:k)='('
        if(text(k:k) == ']' .or. text(k:k) == '}') text(k:k)=')'
      end do

 10   continue
      item = 0
      ilength = len(trim(text))

      if(ilength > 1) then

        do k = 1,(ilength-1)
   
! convert ^ to **, if the user enters ^   MJR unnecessary?
!          if(text(k:k) == '^') then
!            text = text(1:k-1)//'**'//text(k+1:ilength)
!            ilength = ilength + 1
!            item = 1
!          end if

! convert ln to log
! MJR          if(text(k:k+1) == 'ln' .or. text(k:k+1) == 'Ln' .or.text(k:k+1) == 'lN' .or.text(k:k+1) == 'LN') then
! MJR            text = text(1:k-1)//'log'//text(k+2:ilength)
          if(text(k:k+2) == 'ln(' .or. text(k:k+2) == 'Ln(' .or.text(k:k+2) == 'lN(' .or.text(k:k+2) == 'LN(') then
            text = text(1:k-1)//'log('//text(k+3:ilength)
            ilength = ilength + 1
            item = 1
          end if

! convert pi to 3.141592654
! MJR         if(text(k:k+1) == 'pi' .or. text(k:k+1) == 'Pi' .or. text(k:k+1) == 'pI' .or. text(k:k+1) == 'PI') then
! MJR            text = text(1:k-1)//'3.141592654'//text(k+2:ilength)
! MJR           ilength = ilength + 9
! MJR           item = 1
! MJR         end if

! converts comma to full stop, if the user enters comma
          if(text(k:k) == ',') then
            text = text(1:k-1)//'.'//text(k+1:ilength)
            item = 1
          end if

        end do
  
      end if

      if(ilength > 2) then

! penultimate  MJR unnecessary?
!        if(text((ilength-1):(ilength-1)) == '^') then
!          text = text(1:ilength-2)//'**'//text(ilength:ilength)
!          item = 1
!        end if

! penultimate
        if(text((ilength-1):(ilength-1)) == ',') then
          text = text(1:ilength-2)//'.'//text(ilength:ilength)
          item = 1
        end if

      end if

      if(ilength > 1) then
  
! ultimate
        if(text((ilength):(ilength)) == ',') then
          text = text(1:ilength-1)//'.'
          item = 1
        end if
 
      end if

      if(item == 1) goto 10

      return
  
      end subroutine convert_b
!
end module f_e_interpreter
!------------------------------------------------------------------------------
!  1453 lines