                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then ! On top of each other
                               hopp(j,i) = hopp(j,i) + tBA0
                               numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                            else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                               if (Species(i).eq.Species(NList(j,i))) then ! same sublattice
                                  hopp(j,i) = hopp(j,i) + tAA1
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else ! Different sublattice
                                  hopp(j,i) = hopp(j,i) + tAB1
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! BA' (corresponds to g1, not g2)
                               if (Species(i).eq.Species(NList(j,i))) then
                                  hopp(j,i) = hopp(j,i) + tBA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  hopp(j,i) = hopp(j,i) + tBA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aG*2.0_dp*1.1_dp) then ! Second neirest neighbor from central atom in middle of hexagon
                               if (Species(i).eq.Species(NList(j,i))) then
                                  hopp(j,i) = hopp(j,i) + tAA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  hopp(j,i) = hopp(j,i) + tAB2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! BA' (corresponds to g2, not g5)
                               if (Species(i).ne.Species(NList(j,i))) then
                                  hopp(j,i) = hopp(j,i) + tBA5
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  hopp(j,i) = hopp(j,i)
                               end if
                            else
                               hopp(j,i) = hopp(j,i) + 0.0_dp
                            end if
