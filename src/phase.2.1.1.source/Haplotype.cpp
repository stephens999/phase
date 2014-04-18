#include "Haplotype.hpp"

extern int NHAP;

Haplotype::Haplotype(string lt):
  locus_type(lt),
  h(lt.size())
{
  // cout << "constructor called with: " << lt << endl;
  NHAP++;
  //cout << NHAP << endl;
}

Haplotype::Haplotype(Haplotype const & h2):
  h(h2.h),
  locus_type(h2.locus_type)
{
  //  for(int i = 0; i<locus_type.size(); i++)
  // cout << locus_type[i] << endl;

  // cout << "Copy constructor called with: " << locus_type << endl;
  NHAP++;
  //cout << NHAP << endl;
}

Haplotype::Haplotype(const Haplotype & h2, int firstlocus, int lastlocus):
  h(h2.h.begin()+firstlocus,h2.h.begin()+lastlocus),
  locus_type(h2.locus_type.begin()+firstlocus,h2.locus_type.begin()+lastlocus)
{
  NHAP++;
}

// construct a haplotype out of two haplotypes
Haplotype::Haplotype(const Haplotype & h1, const Haplotype & h2):
  h(h1.h),
  locus_type(h1.locus_type+h2.locus_type)
{
  h.insert(h.end(),h2.h.begin(),h2.h.end());
  NHAP++;
}



//  Haplotype::Haplotype(const Haplotype & h2, int firstlocus, int lastlocus):
//    h(h2.h),
//    locus_type(h2.locus_type)
//  {
//    h = vector<int>(h2.h.begin()+firstlocus,h2.h.begin()+lastlocus);
//    locus_type = string(h2.locus_type.begin()+firstlocus,h2.locus_type.begin()+lastlocus);
//  }

Haplotype::~Haplotype()
{
  //cout << "Destructor called" << endl;
  NHAP--;
  //cout << NHAP << endl;
}

Haplotype const & Haplotype::operator=(Haplotype const & h2)
{
  if(this != &h2){
    h = vector<float>(h2.h);
    locus_type = string(h2.locus_type);
  }
  return *this;
}

void Haplotype::print_haplotype(ostream & ostr, const vector<int> * coding) const
{
  for(int locus =0; locus<h.size(); locus++){
    if(locus_type[locus]=='M'){
      if(locus>0)
	if(locus_type[locus-1] == 'S')
	  ostr << ' ';
      ostr << get_allele(locus)-coding[0][locus] << ' ';
    }
    else
      ostr << (char) coding[get_allele(locus)][locus];
  }
}


int NDiff(Haplotype h1, Haplotype h2){ // finds number of diffs
  int d=0;
  for(int l = 0; l < h1.get_nloci(); l++){
    d+= h1.get_allele(l) != h2.get_allele(l); 
  }
  return d;
}
