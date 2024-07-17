/** @type {import('next').NextConfig} */
const nextConfig = {
    rewrites: async () => {
        return [
          {
            source: '/api/:path*',
            destination: 'http://127.0.0.1:8000/api/:path*',
            // destination:
            //   process.env.NODE_ENV === 'development'
            //     ? 'http://127.0.0.1:7000/api/:path*'
            //     : '/api/',
          },
        ]
      },
      // For Docker
      output: "standalone",
};

export default nextConfig;
