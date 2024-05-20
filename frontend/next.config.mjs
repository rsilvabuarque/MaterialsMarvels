/** @type {import('next').NextConfig} */
const nextConfig = {
    rewrites: async () => {
        return [
          {
            source: '/api/:path*',
            destination: 'http://127.0.0.1:7000/api/:path*',
            // destination:
            //   process.env.NODE_ENV === 'development'
            //     ? 'http://127.0.0.1:7000/api/:path*'
            //     : '/api/',
          },
        ]
      },
};

export default nextConfig;
