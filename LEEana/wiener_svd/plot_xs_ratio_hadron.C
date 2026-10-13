TVectorD convert_to_TVectorD(TH1D* h){
	int nbins = h->GetNbinsX();
	TVectorD ret(nbins);
	for (int i=0; i<nbins; i++) {
		ret(i) = h->GetBinContent(i+1);
	}
	return ret;
}

TH1D* convert_to_TH1D(TVectorD v, TString name="", double offset_x=-0.1){
	TH1D* h = new TH1D(name.Data(), "", v.GetNrows(), 0.5 + offset_x, v.GetNrows()+0.5 + offset_x);
	for (int i=0; i<v.GetNrows(); i++) {
		h->SetBinContent(i+1, v(i));
	}
	return h;
}

void put_diag_err(TH1D* h, TMatrixD cov) {
	int nbins = h->GetNbinsX();
	int nele = cov.GetNrows();
	if (nbins != nele or nele==0) {
		cout << "WARNING! Invalid input." << endl;
		return;
	}
	for(int i=0; i<nbins; i++) {
		double err = sqrt(cov(i,i)); 
		h->SetBinError(i+1, err);
	}
}

TMatrixD convert_to_TMatrixD(TH2D* hh){
	int nbinsX = hh->GetNbinsX();
	int nbinsY = hh->GetNbinsY();
	TMatrixD ret(nbinsX, nbinsY); // FIXME: be careful rows and cols
	for (int i=0; i<nbinsX; i++) {
		for (int j=0; j<nbinsY; j++) {
			ret(i,j) = hh->GetBinContent(i+1, j+1);
		}
	}
	return ret;
}



TMatrixD cov_of_Ax(TMatrixD A, TMatrixD covx) { // cov of Ax
	TMatrixD AT = A;
	AT.T(); // invert
	auto ret = A * covx * AT;
	return ret;
}

double sum_of_vector(TVectorD Y){
	double ret=0;
	int nele = Y.GetNrows();
	for (int i=0; i<nele; i++) {
		ret += Y(i);
	}
	return ret;
}

double sum_of_matrix(TMatrixD M){
	double ret=0;
	int nrows = M.GetNrows();
	int ncols = M.GetNcols();
	for (int i=0; i<nrows; i++) {
		for (int j=0; j<ncols; j++) {
			ret += M(i,j);
		}
	}
	return ret;
}

double sum_of_matrix_row(TMatrixD M, int row) {
	double ret=0;
	// int nrows = M.GetNrows();
	int ncols = M.GetNcols();
	for (int j=0; j<ncols; j++) {
		ret += M(row,j);
	}
	return ret;
}

double sum_of_matrix_col(TMatrixD M, int col){
	double ret=0;
	int nrows = M.GetNrows();
	for (int i=0; i<nrows; i++) {
		ret += M(i,col);
	}
	return ret;
}

void convert_to_fractional_uncer(TH1D* h) {
	int nbins = h->GetNbinsX();
	for (int i=0; i<nbins; i++) {
		double content = h->GetBinContent(i+1);
		double error = h->GetBinError(i+1);
		h->SetBinContent(i+1, 0);
		h->SetBinError(i+1, error / content);
	}
}

double GoF(TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD cov){
  TMatrixD md = matrix_pred - matrix_data;
  TMatrixD mdT = matrix_pred - matrix_data;
  mdT.T();
  cov.Invert();
  auto mret = md * cov * mdT;
  cout << mret.GetNrows() << " x " << mret.GetNcols() << " chi2/ndf= " << mret(0,0) << "/" << md.GetNcols() << endl;
  return mret(0,0);
}

void plot_xs_ratio_hadron(){

	// TVectorD X := dsigma/dE
	// TMatrixD covX := covariance of dsigma/dE

	// TMatrixD A = (dE)_i // diagonal matrix
	// TVectorD Y := A*X // dsigma
	// TMatrixD covY :=  A*covX*AT // covariance of dsigma

	// double NT = sum of Y's elements
	// TVectorD Z = (Y_i / NT);
	// TMatrixD covZ = (Mike Shaevitz's equation, drived from covY)

	// TMatrixD B = (1/dE)_i // diagonal matrix
	// TVectorD W = B * Z; // dsigma/dE / sigma -> cancel the normalization uncertainty
	// TMatrixD covW = B * covZ * BT;


	auto file = TFile::Open("output.root");
	auto unf = (TH1D*)file->Get("unf"); // X: = dsigma/dE 
	auto unfcov = (TH2D*)file->Get("unfcov"); // cov(X)

	//
	TVectorD X = convert_to_TVectorD(unf);
	TMatrixD covX = convert_to_TMatrixD(unfcov);

	//
	int nbins = 11;
	double xbins1[] = {0.03, 0.1, 0.15, 0.225, 0.275, 0.336, 0.411, 0.502, 0.614, 0.75, 1.12, 2.5};
	TMatrixD A(nbins, nbins);
	TMatrixD B(nbins, nbins);
	for (int i=0; i<nbins; i++) {
		for (int j=0; j<nbins; j++) {
			if (i==j) {
				A(i,j) = xbins1[i+1] - xbins1[i];
				B(i,j) = 1./A(i,j);
			}
		}
	}

	TVectorD Y = A * X;
	TMatrixD covY = cov_of_Ax(A, covX);

	//
	double NT= sum_of_vector(Y);
	cout << "NT: " << NT << endl;
	int nele = Y.GetNrows();
	TVectorD Z(nele);
	for (int i=0; i<nele; i++) {
		Z(i) = Y(i) / NT;
	} // start to drive covZ with Mike Shaevtiz's formula
	TMatrixD covZ(nele, nele);
	double sum_covY = sum_of_matrix(covY);
	for (int k=0; k<nele; k++) {
		for (int l=0; l<nele; l++) {
			covZ(k,l)  = covY(k,l);
			covZ(k,l) -= Y(l)/NT * sum_of_matrix_row(covY, k);
			covZ(k,l) -= Y(k)/NT * sum_of_matrix_col(covY, l);
			covZ(k,l) += Y(k)*Y(l)/pow(NT,2) * sum_covY;
			covZ(k,l) /= pow(NT,2);
		}
	}


	//
	TVectorD W = B * Z;
	TMatrixD covW = cov_of_Ax(B, covZ);
	// covW.Draw("colz");

	//
	TH1D* hW = convert_to_TH1D(W, "dsigma_dE_sigma");
	put_diag_err(hW, covW);
	hW->SetLineColor(4);
	hW->SetMarkerColor(4);
	hW->GetXaxis()->SetTitle("Energy Bin");
	hW->Draw("e");

	// auto c2 = new TCanvas("c2","",600,450);
	put_diag_err(unf, covX);
	unf->Sumw2();
	unf->Scale(1./NT);
	unf->Draw("e same");

	// convert_to_fractional_uncer(hW);
	// convert_to_fractional_uncer(unf);
	// hW->GetYaxis()->SetRangeUser(-0.5,0.5);
	// hW->Draw("e");
	// unf->Draw("e same");

	//
	double mcbin[11] = {1,2,3,4,5,6,7,8,9,10,11};
	// dsigma/dE/sigma from MC truth
	double mc[11] = {0.941923, 1.3207, 1.39366, 1.52335, 1.74813, 1.73729, 1.34157, 0.888747, 0.523373, 0.237984, 0.050399};
	// dsigma/dE from MC truth
	// double mc[11] = {0.227487, 0.318966, 0.336589, 0.367911, 0.422196, 0.419579, 0.324006, 0.214644, 0.126402, 0.0574763, 0.012172};
	auto gmc = new TGraph(11, mcbin, mc);
	gmc->SetLineStyle(2);
	gmc->SetLineWidth(2);
	gmc->SetLineColor(1);
	gmc->Draw("Lsame");

	TMatrixD m_mc(1,nbins);
	TMatrixD m_data(1,nbins);
	for(int i=0; i<nbins; i++){
	  m_mc(0,i) = mc[i];
	  m_data(0,i) = W(i);
	}
	covW.SetTol(1.e-20); // cannot covW.Invert()?  or we can scale to a larger value, invert and scale back
  	double chi2_mc = GoF(m_mc, m_data, covW);

	// 
	auto lg = new TLegend(0.6,0.6,0.8,0.8);
	// lg->AddEntry(hW, "d#sigma/dE_{#mu}/#sigma")->SetTextColor(4);
	lg->AddEntry(hW, Form("d#sigma/dE_{had}/#sigma, #chi^{2}/ndf: %.1f/%d", chi2_mc,nbins-1) )->SetTextColor(4);
	lg->AddEntry(gmc, "MC Truth", "l");
	// lg->AddEntry(gmc, Form("MC Truth, #chi^{2}/ndf: %.1f/%d",chi2_mc,nbins-1), "l");
	lg->AddEntry(unf,"Unfolded d#sigma/dE_{had} (#times 1/#sigma)");
	lg->Draw();
}
