local gro_read_frame = function(file)
    local title_string = file:read()
    if not title_string then return nil end
    local number_of_atoms = tonumber(file:read())
    if not number_of_atoms then return nil end
    local frame = {}
    frame.natoms = number_of_atoms
    frame.title = title_string
    for k=1,number_of_atoms do
        local l = file:read()
        if not l then return nil end
        frame[#frame+1] = {
            tonumber(l:sub(1,5)),
            l:sub(6,10):gsub("%s+",""),
            l:sub(11,15):gsub("%s+",""),
            tonumber(l:sub(16,20)),
            tonumber(l:sub(21,28)), -- nm
            tonumber(l:sub(29,36)), -- nm
            tonumber(l:sub(37,44)), -- nm
            tonumber(l:sub(45,52)), -- nm/ps
            tonumber(l:sub(53,60)), -- nm/ps
            tonumber(l:sub(61,69)) -- nm/ps
        }
        assert(#frame == k)
    end
    local bx,by,bz = file:read():match('%s*(%S+)%s+(%S+)%s+(%S+)%s*')
    bx,by,bz = tonumber(bx),tonumber(by),tonumber(bz)
    if not (bx and by and bz) then return nil end
    frame.box = {bx,by,bz}
    return frame
end

local split_residues = function(groframe)
    local chain = 0
    local residues = {}
    local id = false
    for k=1,#groframe do
        local resid,resname,atmname,idx,x,y,z = table.unpack(groframe[k])    
        if resid ~= id then
            if resid == 1 then chain = chain + 1 end
            id = resid
            residues[#residues+1] = {resid=resid,resname=resname,chain=chain}
        end
        local res = residues[#residues]
        res[#res+1] = {atmname,x,y,z}
    end
    return residues
end

local coarsemap = {}
coarsemap.BGLC = function(res)
    local ret = {resname=res.resname,resid=res.resid,chain=res.chain}
    local c = {}
    -- map atoms coordinates
    for k=1,#res do
        local atmname,x,y,z,_ = table.unpack(res[k])
        c[atmname] = {x,y,z}
    end
    --for k,v in pairs(c) do print('c',k,v) end
    -- coarsening coordinates
    local set = {{'O2'},{'O3'},{'C6'},{'O6'}}
    for k=1,#set do
        local setk = set[k]
        local n = #setk
        local x,y,z = 0.0,0.0,0.0
        for l=1,n do
            local atm = setk[l]
            local coords = c[atm]
            x = x + coords[1]
            y = y + coords[2]
            z = z + coords[3]
        end
        ret[k] = {'CG'..k,x/n,y/n,z/n}
    end
    return ret
end

local atoms_chain = {}
local coarse = function(residues)
    local c = {}
    local natoms = 0
    local nchains = 0
    local chain = false
    for k=1,#residues do
        local resk = residues[k]
        local resname,chain = resk.resname,resk.chain
        local resc = coarsemap[resname](resk)
        c[#c+1] = resc
        natoms = natoms + #resc
        for l=1,#resc do atoms_chain[#atoms_chain+1] = chain end
    end
    c.natoms = natoms
    print('NATOMS: ', natoms)
    print('NCHAINS:',atoms_chain[#atoms_chain])
    return c
end

local file = assert(io.open('fibril.gro','r'))
local groframe = gro_read_frame(file)
file:close()

local system = coarse(split_residues(groframe))

-----------------------------------------------------------------------------------------------------------------------------------------------------
-- print .gro
-----------------------------------------------------------------------------------------------------------------------------------------------------
local file = assert(io.open('cg-coordinates.gro','w'))
local print = function(str)
    file:write(str,'\n')
end
print('Cellulose coarse-grain chain')
print(system.natoms)
local count = 0
for k=1,#system do
    local res = system[k]
    local resid,resname,chain = res.resid,res.resname,res.chain
    for l=1,#res do
        count = count + 1
        local atmname,x,y,z = table.unpack(res[l])
        res[l][(#(res[l]))+1] = count -- global atom id
        print(string.format('%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f',resid,resname,atmname,count,x,y,z,0.0,0.0,0.0))
    end
end
print('0.000000 0.000000 0.000000')
print''
assert(count == system.natoms)
file:close()

-----------------------------------------------------------------------------------------------------------------------------------------------------
-- print chain topology
-----------------------------------------------------------------------------------------------------------------------------------------------------

local natoms = 4*100

local file = assert(io.open('topol.top','w'))
local print = function(str) file:write(str,'\n') end
--#include "martini_v2.2.itp"
print[[
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
[ moleculetype ]
; name  nrexcl
CG-CELLULOSE    4
]]
print('[ atoms ]')
local map = {CG1='ATX1',CG2='ATX1',CG3='TC1',CG4='ATX1'}
local count = 0
for k=1,#system do
    local res = system[k]
    local resid,resname,chain = res.resid,res.resname,res.chain
    for l=1,#res do
        local atmname,x,y,z,gidx = table.unpack(res[l])
        local martiniatom = map[atmname]
        local charge,mass = 0.0,0.0
		count = count + 1
		if count <= natoms then
    	    print(string.format('%d %s %d %s %s %d',gidx,martiniatom,resid,resname,atmname,gidx))
	    end
	end
end

print''
print('[ bonds ]')
local printbonds = function(a,b,c,d)
    if a <= natoms and b <= natoms and c <= natoms then print(string.format('%d %d 1 %f %f',a,b,c,d)) end
end
for s=0,natoms,4 do
	-- (ALPHA+BETA)/2 GEOMETRY
	-- inter residue
    printbonds(s+1,s+2,0.288,30000.0) -- O2O3 | 2.881358 0.079220 2.532484 3.272505 | 2.882215 0.078287 2.519388 3.304742 | 2.908325 0.080509 2.553345 3.368977
    printbonds(s+1,s+3,0.562,30000.0) -- O2C6 | 5.614538 0.070868 5.234551 5.941634 | 5.621642 0.071645 5.206381 5.976185 | 5.636716 0.070663 5.278480 5.975661
    printbonds(s+2,s+3,0.494,30000.0) -- O3C6 | 4.957965 0.076292 4.555475 5.322707 | 4.958374 0.075792 4.544852 5.335012 | 4.933024 0.078173 4.526765 5.321456 
    printbonds(s+1,s+4,0.642,999999.0) -- O2O6 | 6.502024 0.173141 5.338142 6.969646 | 6.431115 0.212127 5.286103 6.990696 | 6.347560 0.125028 5.356040 6.942769
    printbonds(s+2,s+4,0.576,999999.0) -- O3O6 | 5.611158 0.294130 4.685321 6.507639 | 5.596012 0.328396 4.595749 6.505277 | 6.075645 0.223052 4.630463 6.502710
    printbonds(s+3,s+4,0.143,30000.0) -- C6O6 | 1.430226 0.026871 1.298400 1.560385 | 1.429836 0.026900 1.291097 1.557098 | 1.429587 0.026893 1.283602 1.565965

	-- intra residue
    printbonds(s+1,4+s+1,0.756,30000.0) -- O2+O2 | 7.562013 0.117576 6.427918 8.169572 | 7.567706 0.115508 6.453087 8.167201 | 7.575528 0.120546 6.423901 8.179834
    printbonds(s+2,4+s+2,0.669,30000.0) -- O3+O3 | 6.684076 0.113370 5.731843 7.508379 | 6.699276 0.113105 5.967538 7.577817 | 6.704416 0.123711 5.136170 7.495159
    printbonds(s+3,4+s+3,0.747,30000.0) -- C6+C6 | 7.490481 0.145581 5.699035 8.130223 | 7.480583 0.140999 5.856174 8.189798 | 7.461412 0.128133 5.195501 8.081017

	-- second neighbohrs
    printbonds(s+1,8+s+1,1.044,30000.0) -- O2++O2 | 10.445180 0.196457 8.478990 12.662869 | 10.441422 0.187376 8.885116 12.379725 | 10.422297 0.204867 8.540602 12.222756
    printbonds(s+2,8+s+2,1.044,30000.0) -- O3++O3 | 10.442904 0.155629 9.063269 12.427008 | 10.441671 0.150269 9.376107 12.333647 | 10.425305 0.158044 9.137577 12.196394
    printbonds(s+3,8+s+3,1.044,30000.0) -- C6++C6 | 10.443584 0.210724 7.710785 12.070874 | 10.443072 0.203566 8.579149 11.775724 | 10.423563 0.202031 7.993949 12.046287
end
print''
print[[
#ifdef POSRESCG1CG2
[ position_restraints ]
; atom  type      fx      fy      fz
]]
for k=1,natoms,4 do
    print(string.format('%5d   1  1000  1000  1000',k))
    print(string.format('%5d   1  1000  1000  1000',k+1))
end
print('#endif')
print''

print[[
[ system ]
Title of the system

[ molecules ]
CG-CELLULOSE 36

]]
file:close()












