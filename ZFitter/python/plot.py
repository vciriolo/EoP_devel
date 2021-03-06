import ROOT
from ElectronCategory import cutter,defaultEnergyBranch,ReduceFormula

import errno
import os
def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise

colors = [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kCyan, ROOT.kMagenta, ROOT.kOrange-6, ROOT.kViolet-6, ROOT.kGray+1, ROOT.kOrange-3, ROOT.kViolet-9]
ndraw_entries = 1000000000
#ndraw_entries = 100000

hist_index = 0

def GetTH1(chain, branchname, isMC, binning="", category="", sel="", label="", histname="",
		 usePU=True, useEleIDSF=2, smear=True, scale=True, useR9Weight=False, energyBranchName = None, useLT=True, noWeights = False):
	global hist_index

	if scale and not isMC and CheckForBranch(chain, "scaleEle"):
		branchname = ApplyScaleSmear(branchname, "scaleEle")
	if smear and isMC and CheckForBranch(chain, "smearEle"):
		branchname = ApplyScaleSmear(branchname, "smearEle")

	ActivateBranches(chain, branchname, category, energyBranchName, scale, smear)
	selection, weights = GetSelectionWeights(chain, category, isMC, smear, scale, useR9Weight, usePU, useEleIDSF, useLT, noWeights)
	if(sel!=""):
		selection = ROOT.TCut(str(selection) + "&&" +  sel)
	print "[DEBUG] Getting TH1F of " + branchname + binning + " with " + str(selection) + " and weights" + str(weights)
	if histname:
		chain.Draw(branchname + ">>" + histname, selection * weights, "", ndraw_entries)
	else:
		histname = "hist%d" % hist_index
		chain.Draw(branchname + ">>" + histname + binning, selection * weights, "", ndraw_entries)
		hist_index += 1

	h = ROOT.gDirectory.Get(histname)
	h.SetTitle(label)
	return h

def GetTH1Stack(chain, splits, branchname, isMC, binning="", category="", label="",
		 usePU=True, useEleIDSF=2,smear=False, scale=False, useR9Weight=False, energyBranchName = None, useLT = True):
	''' Unused, use PlotDataMC with stack_mc=True '''
	global hist_index
	ActivateBranches(chain, branchname, category, energyBranchName)
	selection, weights = GetSelectionWeights(chain, category, isMC, smear, scale, useR9Weight, usePU, useEleIDSF, useLT)

	#print "[DEBUG] Getting TH1F of " + branchname + binning + " with " + str(selection) + " and weights" + str(weights)
	
	stack_labels = ["LT_5To75", "LT_75To80", "LT_80To85", "LT_85To90", "LT_90To95", "LT_95To100", "LT_100To200", "LT_200To400", "LT_400To800", "LT_800To2000",]

	nbins, low_b, hi_b = eval(binning)
	# Make histograms
	stack = [ ROOT.TH1F(label + stack_labels[i], stack_labels[i], nbins, low_b, hi_b) for i in range(len(stack_labels))]

	oldfilename = ""
	index = -1
	selector = ROOT.TTreeFormula("selector", selection.GetTitle(), chain)
	weighter = ROOT.TTreeFormula("weighter", weights.GetTitle(), chain)
	brancher = ROOT.TTreeFormula("brancher", branchname, chain)
	maxevents = chain.GetEntries()
	i = 0
	for event in chain:
		if not (i % (maxevents/100)): print 100*i/maxevents, "%"
		i += 1
		selector.UpdateFormulaLeaves()
		sv = selector.EvalInstance()
		if not sv > 0: continue
		filename = event.GetFile().GetName()
		if filename != oldfilename:
			for nlabel, label in enumerate(stack_labels):
				if label in filename:
					index = nlabel
			print index, filename
			oldfilename = filename

		weighter.UpdateFormulaLeaves()
		weight = float(weighter.EvalInstance())

		brancher.UpdateFormulaLeaves()
		value = brancher.EvalInstance()
		stack[index].Fill(value, weight)
		
	return stack

def GetSelectionWeights(chain, category, isMC, smear, scale, useR9Weight, usePU, useEleIDSF, useLT, noWeights):
	selection = ROOT.TCut("")
	if(category):
		try:
			cat1, cat2, common  = category
			selection = cutter.GetCut(cat1 + '-' + common, isMC, 1 , scale) + cutter.GetCut(cat2 + '-' + common, isMC, 2 , scale)
			category = common
		except (TypeError, ValueError):
			selection = cutter.GetCut(category, isMC, 0 , scale and not isMC)

	weights = ROOT.TCut("")
	if(noWeights): return selection, weights
	if(isMC):
		chain.SetBranchStatus("mcGenWeight", 1)
		weights *= "mcGenWeight"
		if(usePU and CheckForBranch(chain, "puWeight")):
			chain.SetBranchStatus("puWeight", 1)
			weights *= "puWeight"
		if(useEleIDSF):
			start = category.find("eleID_") + 6
			end = category.find("-",start)
			if end > 0:
				EleIDSF = "EleIDSF_" + category[start:end]
			else: 
				EleIDSF = "EleIDSF_" + category[start:]

			if CheckForBranch(chain, EleIDSF):
				if useEleIDSF == 2:
					weights *= "{id}[0]*{id}[1]".format(id=EleIDSF)
				elif useEleIDSF == 1:
					weights *= EleIDSF
				chain.SetBranchStatus(EleIDSF, 1)
			else:
				print EleIDSF + " branch not present"


	if useLT and CheckForBranch(chain, "LTweight"):
		chain.SetBranchStatus("LTweight", 1)
		weights *= "LTweight"

	return selection, weights

def ApplyScaleSmear(branchname, scaleOrSmear):
	import re
	branches = ReduceFormula(branchname)
	energy_branches = [b for b in branches if b.startswith("energy")]
	for sub in energy_branches:
		#replace energyBranchname with energyBranchname*scaleEle
		branchname = re.sub(sub + r"(?!\[\d\])", "(" + sub + "*" + scaleOrSmear + ")",branchname)
		#replace energyBranchname[i] with energyBranchname[i]*scaleEle[i]
		branchname = re.sub(sub + r"(\[\d\])", "(" + sub + r"\1" + "*" + scaleOrSmear + r"\1)",branchname)
	invmass_branches = [b for b in branches if b.startswith("invMass")]
	for sub in invmass_branches:
		#replace energyBranchname with energyBranchname*scaleEle
		branchname = re.sub(sub, "(" + sub + "*sqrt({s}[0]*{s}[1]))".format(s=scaleOrSmear), branchname)
	return branchname

def CheckForBranch(chain, branch):
	chain.GetEntry(0)
	try:
		chain.__getattr__(branch)
	except AttributeError:
		return False
	return True


def ActivateBranches(chain, branchnames, category, energyBranchName=None, scale=False, smear=False, debug=False):
	chain.SetBranchStatus("*", 0)
	if energyBranchName:
		cutter.energyBranchName=energyBranchName
	else:
		cutter.energyBranchName=defaultEnergyBranch
		
	branchlist = []
	if category:
		for c in AsList(category):
			branchlist += [str(b) for b in cutter.GetBranchNameNtupleVec(c)]

	branchlist = [x for x in set(branchlist)]
	if scale:
		branchlist += ["scaleEle"]
	if smear:
		branchlist += ["smearEle"]
	for b in branchlist:
		if not CheckForBranch(chain, b): continue
		if debug: print "[Status] Enabling branch:", b
		chain.SetBranchStatus(b, 1)
	for branchname in AsList(branchnames):
		for br in ReduceFormula(branchname):
			if not br: continue
			if not CheckForBranch(chain, br): continue
			if debug: print "[Status] Enabling branch:", br
			chain.SetBranchStatus(br, 1)


def AsList(arg):
	try:
		if 'TH1' in arg.ClassName():
			return [arg,]
		elif "Array" in arg.ClassName():
			return arg
		else:
			print "[WARNING] AsList() encountered unhandled ROOT Class" + arg.ClassName()
			raise
	except:
		if(isinstance(arg,str)):
			return [arg]
		else:
			try:
				iter(arg)
				return arg
			except:
				return [arg]

def ColorMCs(mcs, fill=False):
	mc_hists = AsList(mcs)

	for i, h in enumerate(mc_hists):
		h.SetMarkerStyle(20)
		h.SetMarkerSize(1)
		h.SetMarkerColor(colors[i % len(colors)])
		if fill:
			h.SetFillStyle(1001)
		else:
			h.SetFillStyle(0)
		h.SetFillColorAlpha(colors[i % len(colors)], 1)
		h.SetLineColor(colors[i % len(colors)])
		h.SetLineWidth(2)

def ColorData(data):
	d_hists = AsList(data)

	data_color = [ROOT.kBlack] + colors

	for i, h in enumerate(d_hists):
		h.SetMarkerStyle(20)
		h.SetMarkerSize(1)
		h.SetMarkerColor(data_color[i])

def Normalize(data, mc):
	"""
	Normalize MC histograms to data.Integral()
	if there is no data histogram, normalize to 1
	if there is more than one data histogram
	"""
	data_list = AsList(data)
	mc_list = AsList(mc)
	ndata = len(data_list)
	nmc = len(mc_list)

	if nmc > 0:
                if ndata > 1:
                        for h in data_list:
                                h.Scale(1./h.Integral())
		if ndata:
			for h in mc_list:
				h.Scale(data_list[0].Integral()/h.Integral())
			return
		else:
			for h in mc_list:
				h.Scale(1./h.Integral())
			return
	if ndata > 1 and mc == 0:
		for h in data_list:
			h.Scale(1./h.Integral())
                # for h in mc_list:
		# 	h.Scale(1./h.Integral())
	else:
		print "[WARNING] No normalization defind for (ndata=%d, nmc%d)" % (ndata, nmc)

def PlotDataMC(data, mc, file_path, file_basename, xlabel="", ylabel="",
		ylabel_unit="", logx = False, logy=False, ratio=False, stack_data=False, stack_mc=False,
		x_range = None):

	mc_list = AsList(mc)
	data_list = AsList(data)


	c = ROOT.TCanvas("c", "c", 600, 600)
	c.SetLeftMargin(.14)
	if(ratio):
		c.Divide(1,2)
		c.cd(1)
		ROOT.gPad.SetPad("normal", "normal", 0, .3, 1, 1, ROOT.kWhite, 0, 0);
		ROOT.gPad.SetBottomMargin(0.01)

	nhists = len(data_list) + len(mc_list)
	ncolumns = min(4,nhists)
	nrows = (nhists+1)/ncolumns
	ROOT.gPad.SetTopMargin(.06*nrows)
	ROOT.gPad.SetLogx(logx)
	ROOT.gPad.SetLogy(logy)

	maximum = max( [h.GetMaximum() for h in  data_list] + [h.GetMaximum() for h in mc_list])
	if(logy):
		minimum = 0.001
		maximum *= 5
		ROOT.gPad.SetLogy()
	else:
		minimum = 0.0
		maximum *= 1.2
	
	mcstack = None
	if stack_mc:
		mcstack = ROOT.THStack("mcstack","")
		for h in mc_list:
			mcstack.Add(h)
		mcstack.Draw()

	dstack = None
	if stack_data:
		dstack = ROOT.THStack("dstack","")
		for h in data_list:
			dstack.Add(h)
		dstack.Draw()


	if mcstack:
		h0 = mcstack
	elif dstack:
		h0 = dstack
	elif mc_list:
		h0 = mc_list[0]
	elif data_list:
		h0 = data_list[0]
	else:
		raise Exception("no histograms")

	if not ratio:
		h0.GetXaxis().SetTitle(xlabel)
	else:
		h0.GetXaxis().SetTitleOffset(3)
		h0.GetXaxis().SetLabelOffset(3)

	h0.GetYaxis().SetTitle(ylabel + " /(%.2f %s)" %(h0.GetXaxis().GetBinWidth(2), ylabel_unit))
	h0.GetYaxis().SetRangeUser(minimum, maximum)
	h0.GetYaxis().SetTitleOffset(1.5)

	if x_range:
		h0.GetXaxis().SetRangeUser(*x_range)

	same = ""
	if mcstack:
		mcstack.Draw("hist" + same)
		same = "same"

	#draw stacks
	if dstack:
		dstack.Draw("p" + same)
		same = "same"

	#draw not stacks
	if not stack_mc:
		for h in mc_list:
			h.Draw("hist" + same)
			same = "same" 

	if not stack_data:
		for h in data_list:
			h.Draw("p" + same)
			same = "same" 

	x0 = ROOT.gPad.GetLeftMargin()
	y0 = 1 - ROOT.gPad.GetTopMargin()
	x1 = 1 - ROOT.gPad.GetRightMargin()
	y1 = 1 - ROOT.gPad.GetTopMargin()*.1
	leg = ROOT.TLegend(x0, y0, x1, y1)
	for h in data_list:
		leg.AddEntry(h, h.GetTitle(), "p")
	for h in mc_list:
		leg.AddEntry(h, h.GetTitle(), "l")
	leg.SetNColumns(ncolumns)
	leg.Draw()

	#pv = ROOT.TPaveText(0.25, 0.92, 0.65, .97, "NDC")
	#pv.AddText("CMS Preliminary 2016")
	#pv.SetFillColor(0)
	#pv.SetBorderSize(0)
	#pv.Draw()

	if x_range:
		low,hi = x_range
	else:
		low,hi = h0.GetXaxis().GetXmin(), h0.GetXaxis().GetXmax()
	htemp = ROOT.TH1F("dummy","",1,low, hi)
	if(ratio):
		if not mc_list or not data_list:
			print "[ERROR] Ratio plot with either no MC or no data"
			return
		c.cd(2)
		ROOT.gPad.SetPad("ratio", "ratio", 0.0,  0.0,  1.0,  0.3, ROOT.kWhite, 0, 0 )
		ROOT.gPad.SetTopMargin(0.02)
		ROOT.gPad.SetBottomMargin(0.2)
		ROOT.gPad.SetLogx(logx)
		if(stack_mc):
			ratio_mc = mcstack.GetStack().Last().Clone()
			ratio_mc.SetLineColor(ROOT.kRed)
		else:
			ratio_mc = [h.Clone() for h in mc_list]

		ratio_mc = AsList(ratio_mc)
		htemp.GetXaxis().SetTitle(xlabel)
		htemp.GetYaxis().SetTitle("MC/Data")
		htemp.GetXaxis().SetTitleSize(1.2*.072)
		htemp.GetXaxis().SetLabelSize(.072)
		htemp.GetYaxis().SetTitleSize(1.2*.072)
		htemp.GetYaxis().SetLabelSize(.072)
		try:
			htemp.GetYaxis().SetRangeUser(*ratio)
		except:
			pass
		if x_range:
			htemp.GetXaxis().SetRangeUser(*x_range)
		htemp.SetBinContent(1, 1)
		htemp.SetLineColor(data_list[0].GetMarkerColor())
		htemp.SetLineWidth(2)
		htemp.Draw("hist")
		same = "same"
		for h in ratio_mc:
			h.Divide(data_list[0])
			h.Draw(same)

	mkdir_p(file_path)
	c.SaveAs(file_path + '/' + file_basename + ".png")
	c.SaveAs(file_path + '/' + file_basename + ".pdf")
	c.SaveAs(file_path + '/' + file_basename + ".C")

def NormalizeStack(stack, normalizeTo=1):
	''' Unused replaced with PlotDataMC with stack_mc=True '''

	stack_sum = sum([h.Integral() for h in stack])
	hstack = ROOT.THStack("hs","")
	for h in stack:
		h.Scale(1./stack_sum)
		hstack.Add(h)

	return hstack.GetStack(), hstack

def FoldTH2(h):
	entries = h.GetEntries()
	for i in range(h.GetNbinsX() + 2):
		for j in range(0,i):
			c = h.GetBinContent(i,j)
			co = h.GetBinContent(j,i)
			h.SetBinContent(j,i, c + co)
			h.SetBinContent(i,j, 0)
	h.SetEntries(entries)

def Draw2D(h, label, name, diagBins=None, plotdir="plots/", logx=False, logy=False, logz=False):
	c = ROOT.TCanvas("c","c", 800, 500)

	c.SetRightMargin(.3)

	h.SetStats()
	ROOT.gStyle.SetOptStat(1111111)
	h.SetYTitle("Ele. 1 " + label)
	h.SetXTitle("Ele. 2 " + label)
	h.Draw("colz")

	#print "###", name, h.Integral(0, h.GetNbinsX() + 1, 0, h.GetNbinsY()+1), h.GetEntries()

	c.Modified()
	c.Update()
	ps = h.GetListOfFunctions().FindObject("stats")
	ps.SetX1NDC(0.8);
	c.Modified()
	c.Update()
	if logx: c.SetLogx()
	if logy: c.SetLogy()
	if logz: c.SetLogz()

	if diagBins:
		boxes = []
		for i in xrange(len(diagBins)-1):
			boxes.append(ROOT.TBox( diagBins[i], diagBins[i], diagBins[i+1], diagBins[i+1]))

		for b in boxes:
			b.SetFillStyle(0)
			b.SetLineColor(ROOT.kRed)
			b.Draw()

	c.SaveAs(plotdir + name + ".png")
	c.SaveAs(plotdir + name + ".pdf")
	c.SaveAs(plotdir + name + ".eps")

def makeDiagBins(h, nevents=10000, startbin = None, min_value = None):
	if not startbin:
		xbin = ROOT.Long()
		ybin = ROOT.Long()
		zbin = ROOT.Long()
		maxbin = h.GetMaximumBin()
		h.GetBinXYZ(maxbin, xbin, ybin, zbin)

		startbin = (xbin + ybin)/2
		if h.GetBinContent(startbin, startbin) < h.GetBinContent(startbin + 1, startbin + 1):
			startbin += 1

	if min_value:
		minbin = h.GetXaxis().FindBin(min_value)
	else:
		minbin = 1

	xbins = h.GetXaxis().GetNbins()
	ybins = h.GetYaxis().GetNbins()

	#get center bin
	i = 0
	integral = 0
	while True:
		lowbin = startbin - i 
		hibin = startbin + i 
		if lowbin >= minbin:
			integral += h.Integral(lowbin, lowbin, lowbin, hibin)	# vertical
			did_vert = True
		else:
			lowbin = minbin
			did_vert = False

		if hibin <= xbins:
			#if we did the vertical integral, shift lowbin by 1
			integral += h.Integral(lowbin + did_vert, hibin, hibin, hibin)   # horizontal
		else:
			hibin = xbins

		if integral > nevents: break
		if (hibin == xbins and lowbin == minbin): break
		i += 1

	bins = [lowbin, hibin+1]
	values = [h.GetXaxis().GetBinLowEdge(lowbin), h.GetXaxis().GetBinUpEdge(hibin)]
	integrals = [integral]

	# do top
	i = bins[-1]
	integral = 0
	while i <= xbins and i <= ybins:
		integral += h.Integral(bins[-1], i, i, i)
		#print "top", bins[-1], i, integral,  h.Integral(bins[-1], i, bins[-1], i)
		if integral > nevents or i == xbins:
			bins.append(i+1)
			values.append(h.GetXaxis().GetBinLowEdge(i+1))
			integrals.append(integral)
			integral = 0
		i += 1

	# do bottom
	i = bins[0]-1
	integral = 0
	while i > 0:
		integral += h.Integral(i, i, i, bins[0]-1)
		#print "bot", i, bins[0]-1, integral, h.Integral(i, bins[0]-1, i, bins[0]-1)
		if integral > nevents or i == minbin:
			bins = [i] + bins
			values = [h.GetXaxis().GetBinLowEdge(i)] + values
			integrals = [integral] + integrals
			integral = 0
		i -= 1
	#print values, bins, integrals	
	return values, bins, integrals
