subroutine echo_top_height(r, azi_level0, azi_ref, azi_iter, refl_ref, refl_iter, elev_ref, elev_iter, elev_level0, eth_thld, ethin, eth, nr, na0, nar, nai)
    implicit none
    ! Input variables
    integer(kind=4), intent(in) :: nr, na0, nar, nai
    real(kind=8), intent(in) :: eth_thld, elev_ref, elev_iter, elev_level0
    real(kind=8), intent(in), dimension(nr) :: r
    real(kind=8), intent(in), dimension(na0) :: azi_level0
    real(kind=8), intent(in), dimension(nar) :: azi_ref
    real(kind=8), intent(in), dimension(nai) :: azi_iter
    real(kind=8), intent(in), dimension(nar, nr) :: refl_ref
    real(kind=8), intent(in), dimension(nai, nr) :: refl_iter
    ! Output variables
    real(kind=8), intent(in), dimension(na0, nr) :: ethin
    real(kind=8), intent(out), dimension(na0, nr) :: eth
    
    ! Parameters
    integer(kind=4) :: i, j, nazi_ref, nazi_iter, npos, npos_level0
    real(kind=8), dimension(nr) :: ground_range_reference, ground_range_slice, ground_range_level0
    real(kind=8) :: thetab, thetaa, gr_ref       
    real(kind=8) :: refa, refb, theta_total
    
    eth = ethin
    
    ground_range_reference = r * cosd(elev_ref)
    ground_range_slice = r * cosd(elev_iter)
    ground_range_level0 = r * cosd(elev_level0)
    
    thetab = elev_ref            
    thetaa = elev_iter
    
    do i = 1, na0                   
        nazi_ref = minloc(dabs(azi_level0(i) - azi_ref), 1)
        nazi_iter = minloc(dabs(azi_level0(i) - azi_iter), 1)
        
        do j = 1, nr
            gr_ref = ground_range_reference(j)
            npos = minloc(dabs(ground_range_reference - gr_ref), 1)
            npos_level0 = minloc(dabs(ground_range_level0 - gr_ref), 1)
            
            refb = refl_ref(nazi_ref, j)
            refa = refl_iter(nazi_iter, npos)
            if(ISNAN(refa)) then
                refa = -2
            endif
            
            if ((npos >= nr).AND.(dabs(ground_range_slice(npos) - gr_ref) < 500)) then
                if (refb >= eth_thld) then
                    eth(i, npos_level0) = thetab + 0.5
                    continue
                endif
            endif
            
            if((refb > eth_thld).AND.(refa <= eth_thld)) then
                if (thetaa < 90) then
                    theta_total = (eth_thld - refa) * (thetab - thetaa) / (refb - refa) + thetab
                else
                    theta_total = thetab + 0.5
                endif
                
                eth(i, npos_level0) = theta_total                
            endif
        enddo
    enddo                
    
    return
end subroutine echo_top_height


subroutine grid_radar(cloudtop, eth_out, xgrid, ygrid, xradar, yradar, theta3db, nx, ny, nxrad, nyrad)
    implicit none
    ! Input variables
    integer(kind=4), intent(in) :: nx, ny, nxrad, nyrad
    real(kind=8), intent(in) :: theta3db
    real(kind=8), intent(in), dimension(nx) :: xgrid
    real(kind=8), intent(in), dimension(ny) :: ygrid
    real(kind=8), intent(in), dimension(nyrad, nxrad) :: xradar
    real(kind=8), intent(in), dimension(nyrad, nxrad) :: yradar
    real(kind=8), intent(in), dimension(nyrad, nxrad) :: cloudtop
    ! Output variables
    real(kind=8), intent(out), dimension(ny, nx) :: eth_out
    
    ! Parameters
    integer(kind=4) :: cnt
    integer(kind=4) :: i, j, k, l
    integer(kind=4) :: stx, sty, edx, edy
    real(kind=8) :: xi, yi, xr, yr
    real(kind=8) :: width, rmax
    real(kind=8) :: zmax
    
    rmax = maxval(xradar)
    
    do i = 1, nx
        do j = 1, ny
            cnt = 0
            zmax = 0.0
            xi = xgrid(i)
            yi = ygrid(j)                        
            
            if (xi ** 2 + yi ** 2 > rmax ** 2) then
                continue
            endif
            
            width = 0.5 * (sqrt(xi ** 2 + yi ** 2) * theta3db * 3.14 / 180.0)            
      
            do k = 1, nxrad
                do l = 1, nyrad
                    xr = xradar(l, k)
                    yr = yradar(l, k)
                    
                    if ((xr >= xi - width).AND.(xr < xi + width).AND.(yr >= yi - width).AND.(yr < yi + width)) then
                        if (cloudtop(l, k) > 0) then
                            zmax = zmax + cloudtop(l, k)
                            cnt = cnt + 1
                        endif
                        continue
                    endif
                enddo
            enddo
            
            eth_out(j, i) = zmax / cnt
        enddo
    enddo
    
    return
    
end subroutine grid_radar