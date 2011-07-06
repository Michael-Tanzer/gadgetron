#ifndef CUCGPRECONDITIONER_H
#define CUCGPRECONDITIONER_H
#pragma once

#include "cuNDArray.h"

template <class T> class cuCGPreconditioner
{
 public:
  cuCGPreconditioner( int device = -1 )
  {
    if( device<0 ){
      if( cudaGetDevice( &device_ ) != cudaSuccess ){
	std::cerr << "cuCGPreconditioner: unable to get current device." << std::endl ;
	device_ = 0;
      }      
    }
    else
      device_ = device;    
  }
  
  virtual ~cuCGPreconditioner() {}
  virtual int apply(cuNDArray<T>* in, cuNDArray<T>* out) = 0;

  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }

protected:
  virtual int set_device()
  {
    if( cudaGetDevice( &old_device_ ) != cudaSuccess ){
      std::cerr << "cuCGPreconditioner::set_device: unable to get current device." << std::endl ;
      return -1;
    }
    if( device_ != old_device_ && cudaSetDevice(device_) != cudaSuccess) {
      std::cerr << "cuCGPreconditioner::set_device: unable to set device " << device_ << std::endl;
      return -1;
    }
    return 0;
  }

  virtual int restore_device()
  {
    if( device_ != old_device_ && cudaSetDevice(old_device_) != cudaSuccess) {
      std::cerr << "cuCGMatrixOperator::restore_device: unable to set device " << old_device_ << std::endl;
      return -1;
    }
    return 0;
  }

 protected:
  int device_;
  
 private:
  int old_device_;

};

#endif //CUCGPRECONDITIONER_H
