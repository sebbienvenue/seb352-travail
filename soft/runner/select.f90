! numerical recipes chapter 8.5
      function select1(k,n,arr)

      integer k,n
      real*8 select1,arr(n)
      integer i,ir,j,l,mid
      real*8 a,temp

      l=1
      ir=n
 1    if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
        endif
        select1=arr(k)
        return
      else
!!        write(*,*)'else ',ir,l
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
 3      continue
          i=i+1
        if(arr(i).lt.a) goto 3
 4      continue
          j=j-1
        if(arr(j).gt.a) goto 4
        if(j.lt.i) goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
 5      arr(l+1)=arr(j)
        arr(j)=a
        if(j.ge.k)ir=j-1
        if(j.le.k)l=i
      endif
      goto 1
      end











