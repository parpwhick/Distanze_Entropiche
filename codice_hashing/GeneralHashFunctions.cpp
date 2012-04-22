
inline unsigned int RSHash(const unsigned char* str, std::size_t len)
{
   unsigned int b    = 378551;
   unsigned int a    = 63689;
   unsigned int hash = 0;

   for(std::size_t i = 0; i < len; i++)
   {
      hash = hash * a + str[i];
      a    = a * b;
   }

   return hash;
}
/* End Of RS Hash Function */


inline unsigned int JSHash(const unsigned char* str, std::size_t len)
{
   unsigned int hash = 1315423911;

   for(std::size_t i = 0; i < len; i++)
   {
      hash ^= ((hash << 5) + str[i] + (hash >> 2));
   }

   return hash;
}
/* End Of JS Hash Function */


inline unsigned int PJWHash(const unsigned char* str, std::size_t len)
{
   unsigned int BitsInUnsignedInt = (unsigned int)(sizeof(unsigned int) * 8);
   unsigned int ThreeQuarters     = (unsigned int)((BitsInUnsignedInt  * 3) / 4);
   unsigned int OneEighth         = (unsigned int)(BitsInUnsignedInt / 8);
   unsigned int HighBits          = (unsigned int)(0xFFFFFFFF) << (BitsInUnsignedInt - OneEighth);
   unsigned int hash              = 0;
   unsigned int test              = 0;

   for(std::size_t i = 0; i < len; i++)
   {
      hash = (hash << OneEighth) + str[i];

      if((test = hash & HighBits)  != 0)
      {
         hash = (( hash ^ (test >> ThreeQuarters)) & (~HighBits));
      }
   }

   return hash;
}
/* End Of  P. J. Weinberger Hash Function */


inline unsigned int ELFHash(const unsigned char* str, std::size_t len)
{
   unsigned int hash = 0;
   unsigned int x    = 0;

   for(std::size_t i = 0; i < len; i++)
   {
      hash = (hash << 4) + str[i];
      if((x = hash & 0xF0000000L) != 0)
      {
         hash ^= (x >> 24);
      }
      hash &= ~x;
   }

   return hash;
}
/* End Of ELF Hash Function */


inline unsigned int BKDRHash(const unsigned char* str, std::size_t len)
{
   unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
   unsigned int hash = 0;

   for(std::size_t i = 0; i < len; i++)
   {
      hash = (hash * seed) + str[i];
   }

   return hash;
}
/* End Of BKDR Hash Function */


inline unsigned int SDBMHash(const unsigned char* str, std::size_t len)
{
   unsigned int hash = 0;

   for(std::size_t i = 0; i < len; i++)
   {
      hash = str[i] + (hash << 6) + (hash << 16) - hash;
   }

   return hash;
}
/* End Of SDBM Hash Function */


inline unsigned int DJBHash(const unsigned char* str, std::size_t len)
{
   unsigned int hash = 5381;

   for(std::size_t i = 0; i < len; i++)
   {
      hash = ((hash << 5) + hash) + str[i];
   }

   return hash;
}
/* End Of DJB Hash Function */


inline unsigned int DEKHash(const unsigned char* str, std::size_t len)
{
   unsigned int hash = static_cast<unsigned int>(len);

   for(std::size_t i = 0; i < len; i++)
   {
      hash = ((hash << 5) ^ (hash >> 27)) ^ str[i];
   }

   return hash;
}
/* End Of DEK Hash Function */


inline unsigned int BPHash(const unsigned char* str, std::size_t len)
{
   unsigned int hash = 0;
   for(std::size_t i = 0; i < len; i++)
   {
      hash = hash << 7 ^ str[i];
   }

   return hash;
}
/* End Of BP Hash Function */


inline unsigned int FNVHash(const unsigned char* str, std::size_t len)
{
   const unsigned int fnv_prime = 0x811C9DC5;
   unsigned int hash = 0;
   for(std::size_t i = 0; i < len; i++)
   {
      hash *= fnv_prime;
      hash ^= str[i];
   }

   return hash;
}
/* End Of FNV Hash Function */


inline unsigned int APHash(const unsigned char* str, std::size_t len)
{
   unsigned int hash = 0xAAAAAAAA;

   for(std::size_t i = 0; i < len; i++)
   {
      hash ^= ((i & 1) == 0) ? (  (hash <<  7) ^ str[i] * (hash >> 3)) :
                               (~((hash << 11) + (str[i] ^ (hash >> 5))));
   }

   return hash;
}
/* End Of AP Hash Function */
