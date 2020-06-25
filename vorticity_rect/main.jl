using Plots



function gauss_seidel(psi, lat, h)
   ERROR = 2
   while ERROR > 0.001
        ERROR = 0
        for i = 2:size(psi,1) - 1
            for j = 2:size(psi,2) - 1
                old_psi = psi[i,j]
                psi[i,j] = 1/4*(psi[i-1,j].+psi[i+1,j]+psi[i,j-1]+psi[i,j+1] + h^2*lat[i,j])
                ERROR += abs(old_psi-psi[i,j])
            end
        end
   plt = contourf(psi)
   display(plt)
   end
end

function main()
   lat_size = 25
   lat = zeros((lat_size,lat_size))
   #lat[5:10,5:10] .= -0.5
   #lat[15:20,15:20] .= -0.5
    
   lat[15:20,15:20] .= 0.25
   lat[5:10,15:20] .= -0.5
   
   lat[5:10,5:10] .= 0.5
   lat[1:6,5:10] .= -1.25
    



   psi = ones((lat_size,lat_size))*0.5
   psi[:,1] .= 0.0
   psi[:,25] .= 0.0
   psi[1,:] .= 0.0
   psi[25,:] .= 0.0
   h = 0.25
   gr() 
   gauss_seidel(psi,lat,h)

end

main()
