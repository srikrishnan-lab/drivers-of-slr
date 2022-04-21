using Mimi
using MimiDICE2013
using MimiSNEASY
#using MimiBRICK

DICE = MimiDICE2013.get_model()
run(DICE)

tot_emissions = getdataframe(DICE, :emissions=>(:CCA, :E, :EIND))
#= 
- CCA is Cumulative indiustrial emissions
- E is Total CO2 emissions (GtCO2 per year)
- EIND is Industrial emissions (GtCO2 per year) 
=#
print(tot_emissions)

co2_emissions = getdataframe(DICE, :emissions=>:E)[60, :]
print(co2_emissions) # just prints total co2 emissions from last year (2305)

# should probably be SNEASY-BRICK
SNEASY = MimiSNEASY.get_model()
run(SNEASY)

#= notes:
- issues w remote access workstation
- how to get BRICK? (issue w private repo?)
- doing Julia env correctly?
    -- is there a way to keep the julia module loaded?
- stuck on generating random parameters
- currently: how to get emissions input to SNEASY-BRICK
=#
