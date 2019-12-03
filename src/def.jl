# define the type of variables
export
  costDataType,GenericData,ParsedData

type costDataType
  model :: Int64
  upCost :: Float64
  downCost :: Float64
  n :: Int64
  params :: Array{Float64,1}
end

type uncertainData
  GPmax :: Float64
  GQmax :: Float64
  GPmin :: Float64
  GQmin :: Float64
  GP0 :: Float64
  GQ0 :: Float64

  WPmax :: Float64
  WQmax :: Float64

  RESPmax :: Float64
  RESQmax :: Float64
  RESPmin :: Float64
  RESQmin :: Float64
  RESP0 :: Float64
  RESQ0 :: Float64

  rPmax :: Float64
  rPmin :: Float64
  rQmax :: Float64
  rQmin :: Float64
end

type fixedData
  baseMVA :: Float64
  bType :: Dict{Int64,Any}

  IDList :: Array{Int64,1}
  genIDList :: Array{Int64,1}
  brList :: Array{Any,1}
  brRev :: Dict{Any,Any}

  Loc :: Dict{Int64,Any}
  LocRev :: Dict{Int64,Any}
  Vmax :: Dict{Int64,Any}
  Vmin :: Dict{Int64,Any}
  Pmax :: Dict{Int64,Any}
  Pmin :: Dict{Int64,Any}
  Qmax :: Dict{Int64,Any}
  Qmin :: Dict{Int64,Any}
  gs :: Dict{Int64,Any}
  bs :: Dict{Int64,Any}
  Vmag :: Dict{Int64,Any}
  Vang :: Dict{Int64,Any}
  Pd :: Dict{Int64,Any}
  Qd :: Dict{Int64,Any}
  Pg :: Dict{Int64,Any}
  Qg :: Dict{Int64,Any}

  g :: Dict{Tuple{Int64,Int64,Int64},Any}
  b :: Dict{Tuple{Int64,Int64,Int64},Any}
  bc :: Dict{Tuple{Int64,Int64,Int64},Any}
  θmax :: Dict{Tuple{Int64,Int64,Int64},Any}
  θmin :: Dict{Tuple{Int64,Int64,Int64},Any}
  rateA :: Dict{Tuple{Int64,Int64,Int64},Any}
  τ1 :: Dict{Tuple{Int64,Int64,Int64},Any}
  τ2 :: Dict{Tuple{Int64,Int64,Int64},Any}
  σ :: Dict{Tuple{Int64,Int64,Int64},Any}

  cp :: Dict{Any,Any}
  cq :: Dict{Any,Any}

  busInd :: Dict{Any,Any}
end
